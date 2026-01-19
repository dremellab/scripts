#!/usr/bin/env python3
"""
download_sra_from_gse.py

Download SRA runs for a GEO Series (GSE...); optionally submit one SLURM job per run.

Local mode (default):
  prefetch -> fasterq-dump -> pigz (default ON)

SLURM mode (--slurm):
  Writes one sbatch script per run and submits it via sbatch.
  Each job runs:
    prefetch -> fasterq-dump -> pigz (default ON)
  and (optionally) cleans up .sra and empty SRR folders.

Metadata:
  --write-metadata:
    Writes these into --outdir:
      <GSE>_series_matrix.txt
      <GSE>_runinfo.csv
      <GSE>_merged_metadata.tsv   (one row per Run)
    IMPORTANT FIX: RunInfo 'Sample' is often SRS..., not GSM....
    We find GSM by scanning each RunInfo row for a token matching GSM\\d+.

  --metadata-only:
    Only write metadata and exit (no downloads and no SLURM submission).

SLURM defaults requested:
  --slurm-partition  standard
  --slurm-account    dremel_lab
  --slurm-time       12:00:00
  --slurm-mem        48G
  --threads          8 (also cpus-per-task)
  job-name           SRR accession
  no qos option
  each sbatch script activates mamba env sra-tools by default

Example:
  # metadata only
  python download_sra_from_gse.py GSE59717 -o /path/to/GSE59717 --metadata-only

  # local download + metadata
  python download_sra_from_gse.py GSE59717 -o /path/to/GSE59717 --write-metadata

  # submit one job per SRR + metadata
  python download_sra_from_gse.py GSE59717 -o /path/to/GSE59717 --write-metadata --slurm
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import time
from typing import Dict, List, Optional, Set, Tuple

import requests

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


# -------------------- EUTILS HELPERS --------------------

def eutils_get(path: str, params: dict, timeout: int = 60) -> str:
    r = requests.get(f"{EUTILS_BASE}/{path}", params=params, timeout=timeout)
    r.raise_for_status()
    return r.text


def esearch(db: str, term: str, retmax: int = 50) -> List[str]:
    r = requests.get(
        f"{EUTILS_BASE}/esearch.fcgi",
        params={"db": db, "term": term, "retmax": retmax, "retmode": "json"},
        timeout=60,
    )
    r.raise_for_status()
    return r.json().get("esearchresult", {}).get("idlist", [])


def elink(dbfrom: str, db: str, ids: List[str]) -> List[str]:
    if not ids:
        return []
    xml = eutils_get(
        "elink.fcgi",
        {"dbfrom": dbfrom, "db": db, "id": ",".join(ids), "retmode": "xml"},
    )
    return re.findall(r"<Id>(\d+)</Id>", xml)


def efetch_runinfo(sra_uids: List[str]) -> str:
    if not sra_uids:
        return ""
    return eutils_get(
        "efetch.fcgi",
        {"db": "sra", "id": ",".join(sra_uids), "rettype": "runinfo", "retmode": "text"},
        timeout=120,
    )


# -------------------- PARSING --------------------

def split_csv(line: str) -> List[str]:
    """Minimal CSV split (handles quotes + doubled quotes)."""
    out: List[str] = []
    cur: List[str] = []
    quoted = False
    i = 0
    while i < len(line):
        c = line[i]
        if c == '"':
            if quoted and i + 1 < len(line) and line[i + 1] == '"':
                cur.append('"')
                i += 1
            else:
                quoted = not quoted
        elif c == "," and not quoted:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(c)
        i += 1
    out.append("".join(cur))
    return out


def parse_runs(runinfo: str) -> Set[str]:
    runs: Set[str] = set()
    if not runinfo.strip():
        return runs

    lines = runinfo.splitlines()
    if not lines:
        return runs

    header = lines[0].split(",")

    try:
        idx = header.index("Run")
        for line in lines[1:]:
            cols = split_csv(line)
            if idx < len(cols) and re.fullmatch(r"[SED]RR\d+", cols[idx]):
                runs.add(cols[idx])
    except ValueError:
        for line in lines[1:]:
            runs.update(re.findall(r"\b[SED]RR\d+\b", line))

    return runs


def strip_wrapping_quotes(s: str) -> str:
    s = s.strip()
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        return s[1:-1]
    return s


# -------------------- UTILS --------------------

def check_tool(name: str) -> None:
    if not shutil.which(name):
        raise RuntimeError(f"Required tool not found on PATH: {name}")


def run_cmd(cmd: List[str]) -> None:
    print("[cmd]", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def sh_quote(s: str) -> str:
    """Minimal safe shell quoting for paths/args."""
    return "'" + s.replace("'", "'\"'\"'") + "'"


def cleanup_empty_run_dir(outdir: str, run_acc: str) -> None:
    """Remove an empty <outdir>/<run_acc>/ directory if it exists."""
    run_dir = os.path.join(outdir, run_acc)
    try:
        if os.path.isdir(run_dir) and not os.listdir(run_dir):
            os.rmdir(run_dir)
    except OSError:
        pass


# -------------------- CORE LOGIC --------------------

def gse_to_sra_uids(gse: str) -> List[str]:
    geo_uids = esearch("gds", f"{gse}[Accession]")
    if not geo_uids:
        raise RuntimeError(f"No GEO record found for {gse}")

    sra_uids = elink("gds", "sra", geo_uids)
    if not sra_uids:
        raise RuntimeError(
            f"No SRA links found for {gse}. (Sometimes GEO links are indirect or absent.)"
        )
    return sra_uids


def gse_to_runs(gse: str) -> List[str]:
    sra_uids = gse_to_sra_uids(gse)
    runinfo = efetch_runinfo(sra_uids)
    runs = sorted(parse_runs(runinfo))
    if not runs:
        raise RuntimeError(f"Failed to extract SRR/ERR/DRR runs for {gse}")
    return runs


# -------------------- METADATA (Series Matrix + RunInfo + Merge) --------------------

def geo_matrix_url(gse: str) -> str:
    m = re.fullmatch(r"GSE(\d+)", gse)
    if not m:
        raise ValueError(f"Not a GSE accession: {gse}")
    digits = m.group(1)
    family = f"GSE{digits[:-3]}nnn"  # e.g. 59717 -> GSE597nnn
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{family}/{gse}/matrix/{gse}_series_matrix.txt.gz"


def download_series_matrix(gse: str, outdir: str) -> str:
    os.makedirs(outdir, exist_ok=True)
    url = geo_matrix_url(gse)
    out_path = os.path.join(outdir, f"{gse}_series_matrix.txt")

    print(f"[info] Downloading series matrix: {url}")
    r = requests.get(url, stream=True, timeout=120)
    r.raise_for_status()

    with gzip.GzipFile(fileobj=r.raw) as gz_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(gz_in, f_out)

    print(f"[info] Wrote: {out_path}")
    return out_path


def write_runinfo_csv(gse: str, outdir: str, sra_uids: List[str]) -> str:
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, f"{gse}_runinfo.csv")
    runinfo = efetch_runinfo(sra_uids)
    if not runinfo.strip():
        raise RuntimeError(f"Empty runinfo returned for {gse} (SRA UIDs: {len(sra_uids)})")

    with open(out_path, "w", encoding="utf-8", newline="") as f:
        f.write(runinfo)

    print(f"[info] Wrote: {out_path}")
    return out_path


def parse_series_matrix(matrix_path: str) -> Dict[str, Dict[str, str]]:
    """Return gsm -> dict of normalized GEO sample_* fields."""
    gsm_ids: List[str] = []
    raw_fields: Dict[str, List[str]] = {}

    with open(matrix_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.startswith("!Sample_"):
                continue
            parts = line.split("\t")
            key = parts[0]
            vals = [strip_wrapping_quotes(v) for v in parts[1:]]

            if key == "!Sample_geo_accession":
                gsm_ids = vals
            raw_fields[key] = vals

    if not gsm_ids:
        raise RuntimeError(f"Could not find !Sample_geo_accession in {matrix_path}")

    per_gsm: Dict[str, Dict[str, str]] = {gsm: {} for gsm in gsm_ids}

    for key, vals in raw_fields.items():
        if not key.startswith("!Sample_"):
            continue
        suffix = key[len("!Sample_"):]  # e.g. title
        norm_key = "sample_" + suffix.lower()  # e.g. sample_title
        for i, gsm in enumerate(gsm_ids):
            per_gsm[gsm][norm_key] = vals[i] if i < len(vals) else ""

    return per_gsm


def parse_runinfo_csv(runinfo_path: str) -> List[Dict[str, str]]:
    with open(runinfo_path, "r", encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f))


def find_gsm_in_runinfo_row(row: Dict[str, str]) -> Tuple[str, str]:
    """Return (GSM accession, column name where found). If none, ('','').

    We scan ALL values for a token matching GSM\\d+.
    """
    for col, v in row.items():
        if not v:
            continue
        m = re.search(r"\bGSM\d+\b", str(v))
        if m:
            return m.group(0), col
    return "", ""


def write_merged_metadata(gse: str, outdir: str, matrix_path: str, runinfo_path: str) -> str:
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, f"{gse}_merged_metadata.tsv")

    per_gsm = parse_series_matrix(matrix_path)
    run_rows = parse_runinfo_csv(runinfo_path)

    sample_keys = sorted({k for d in per_gsm.values() for k in d.keys()})

    run_cols = [
        "Run",
        "Sample",
        "BioSample",
        "Experiment",
        "SRAStudy",
        "LibraryLayout",
        "LibraryStrategy",
        "LibrarySource",
        "LibrarySelection",
        "Instrument",
        "avgLength",
        "bases",
        "size_MB",
    ]

    header = run_cols + ["geo_gsm", "geo_gsm_source"] + sample_keys

    missing = 0
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        for rr in run_rows:
            gsm, src = find_gsm_in_runinfo_row(rr)
            geo = per_gsm.get(gsm, {}) if gsm else {}
            if not gsm:
                missing += 1

            out_row: Dict[str, str] = {}
            for c in run_cols:
                out_row[c] = (rr.get(c) or "").strip()

            out_row["geo_gsm"] = gsm
            out_row["geo_gsm_source"] = src

            for k in sample_keys:
                out_row[k] = geo.get(k, "")

            w.writerow(out_row)

    if missing:
        print(f"[warn] {missing} run(s) had no GSM found in RunInfo; GEO sample_* fields will be blank for those rows.")
    print(f"[info] Wrote: {out_path}")
    return out_path


def write_metadata_bundle(gse: str, outdir: str) -> Tuple[str, str, str]:
    sra_uids = gse_to_sra_uids(gse)
    matrix_path = download_series_matrix(gse, outdir)
    runinfo_path = write_runinfo_csv(gse, outdir, sra_uids)
    merged_path = write_merged_metadata(gse, outdir, matrix_path, runinfo_path)
    return matrix_path, runinfo_path, merged_path


# -------------------- DOWNLOAD / CONVERT --------------------

def compress_fastqs(outdir: str, run_acc: str, threads: int) -> None:
    fastqs = [
        f for f in os.listdir(outdir)
        if f.startswith(run_acc) and f.endswith(".fastq")
    ]
    if not fastqs:
        print(f"[warn] No FASTQs found for {run_acc} to compress")
        return
    run_cmd(["pigz", "-p", str(int(threads))] + [os.path.join(outdir, f) for f in fastqs])


def cleanup_sra(outdir: str, run_acc: str) -> None:
    candidates = [
        os.path.join(outdir, f"{run_acc}.sra"),
        os.path.join(outdir, run_acc, f"{run_acc}.sra"),
    ]
    for p in candidates:
        if os.path.exists(p):
            try:
                os.remove(p)
            except OSError:
                pass


def process_run_local(
    *,
    run_acc: str,
    outdir: str,
    threads: int,
    keep_sra: bool,
    pigz_on: bool,
    pigz_threads: int,
    prefetch_extra: List[str],
    fasterq_extra: List[str],
) -> None:
    run_cmd(["prefetch", run_acc, "-O", outdir] + prefetch_extra)

    run_cmd(
        [
            "fasterq-dump",
            run_acc,
            "--outdir", outdir,
            "--threads", str(int(threads)),
            "--split-files",
        ] + fasterq_extra
    )

    if pigz_on:
        compress_fastqs(outdir, run_acc, pigz_threads)

    if not keep_sra:
        cleanup_sra(outdir, run_acc)

    # remove empty SRR directory if prefetch left it behind
    cleanup_empty_run_dir(outdir, run_acc)


# -------------------- SLURM SUBMISSION --------------------

def write_sbatch_script(
    *,
    run_acc: str,
    outdir: str,
    logdir: str,
    cpus: int,
    mem: str,
    time_limit: str,
    partition: str,
    account: str,
    keep_sra: bool,
    pigz_on: bool,
    pigz_threads: int,
    prefetch_extra: List[str],
    fasterq_extra: List[str],
    prolog: str,
) -> str:
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)

    sb_lines = [
        "#!/usr/bin/env bash",
        "#SBATCH --export=NONE",
        f"#SBATCH --job-name={run_acc}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --account={account}",
        f"#SBATCH --cpus-per-task={int(cpus)}",
        f"#SBATCH --mem={mem}",
        f"#SBATCH --time={time_limit}",
        f"#SBATCH --output={os.path.join(logdir, run_acc)}.%j.out",
        f"#SBATCH --error={os.path.join(logdir, run_acc)}.%j.err",
        "",
    ]

    prefetch_cmd = " ".join(
        ["prefetch", run_acc, "-O", sh_quote(outdir)] + [sh_quote(x) for x in prefetch_extra]
    )
    fasterq_cmd = " ".join(
        [
            "fasterq-dump",
            run_acc,
            "--outdir", sh_quote(outdir),
            "--threads", str(int(cpus)),
            "--split-files",
        ] + [sh_quote(x) for x in fasterq_extra]
    )

    if pigz_on:
        pigz_block = "\n".join(
            [
                "shopt -s nullglob",
                f"files=({sh_quote(outdir)}/{run_acc}*.fastq)",
                "if [ ${#files[@]} -gt 0 ]; then",
                f"  pigz -p {int(pigz_threads)} \"${{files[@]}}\"",
                "else",
                f"  echo \"[warn] No FASTQs found to pigz for {run_acc}\"",
                "fi",
            ]
        )
    else:
        pigz_block = 'echo "[info] pigz disabled"'

    if keep_sra:
        cleanup_block = 'echo "[info] keeping .sra"'
    else:
        cleanup_block = "\n".join(
            [
                f"rm -f {sh_quote(os.path.join(outdir, run_acc + '.sra'))} || true",
                f"rm -f {sh_quote(os.path.join(outdir, run_acc, run_acc + '.sra'))} || true",
            ]
        )

    # Remove empty run folder after cleanup
    empty_dir_cleanup = "\n".join(
        [
            f"run_dir={sh_quote(os.path.join(outdir, run_acc))}",
            "if [ -d \"$run_dir\" ]; then",
            "  if [ -z \"$(ls -A \"$run_dir\" 2>/dev/null)\" ]; then",
            "    rmdir \"$run_dir\" || true",
            "  fi",
            "fi",
        ]
    )

    body = "\n".join(
        [
            "set -euo pipefail",
            "",
            "# Prolog: activate tools/env on compute nodes",
            prolog.strip() if prolog.strip() else "true",
            "",
            'echo "[info] Host: $(hostname)"',
            'echo "[info] Start: $(date)"',
            f'echo "[info] Run: {run_acc}"',
            f"mkdir -p {sh_quote(outdir)}",
            "",
            f'echo "[info] {prefetch_cmd}"',
            prefetch_cmd,
            "",
            f'echo "[info] {fasterq_cmd}"',
            fasterq_cmd,
            "",
            pigz_block,
            "",
            cleanup_block,
            "",
            empty_dir_cleanup,
            "",
            'echo "[info] Done: $(date)"',
        ]
    )

    return "\n".join(sb_lines) + body + "\n"


def submit_srr_as_slurm_job(*, run_acc: str, jobs_dir: str, sbatch_script: str, dry_run: bool) -> Optional[str]:
    os.makedirs(jobs_dir, exist_ok=True)
    script_path = os.path.join(jobs_dir, f"{run_acc}.sbatch")
    with open(script_path, "w", encoding="utf-8") as f:
        f.write(sbatch_script)

    if dry_run:
        print(f"[slurm] (dry-run) wrote: {script_path}")
        return None

    print(f"[slurm] submitting: {script_path}")
    p = subprocess.run(["sbatch", script_path], capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            f"sbatch failed for {run_acc}:\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )

    m = re.search(r"Submitted batch job\s+(\d+)", p.stdout)
    job_id = m.group(1) if m else None
    if job_id:
        print(f"[slurm] {run_acc} -> job {job_id}")
    else:
        print(f"[slurm] {run_acc} submitted; could not parse job id. sbatch said: {p.stdout.strip()}")
    return job_id


# -------------------- MAIN --------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="Download SRA runs for a GSE; optionally submit one SLURM job per run.")

    ap.add_argument("gse", help="GEO Series accession, e.g. GSE59717")
    ap.add_argument("-o", "--outdir", default="sra_downloads", help="Output directory")

    ap.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Threads for fasterq-dump (and default pigz threads). Default: 8",
    )

    ap.add_argument("--limit", type=int, default=0, help="Only process first N runs (0 = all)")
    ap.add_argument("--keep-sra", action="store_true", help="Keep .sra after conversion")
    ap.add_argument("--sleep", type=float, default=0.0, help="Sleep between runs (local mode only)")

    # pigz options (default ON)
    ap.add_argument("--no-pigz", dest="pigz", action="store_false", help="Disable pigz compression")
    ap.add_argument("--pigz-threads", type=int, default=None, help="pigz threads (default: --threads)")
    ap.set_defaults(pigz=True)

    ap.add_argument("--prefetch-extra", default="", help='Extra args for prefetch, e.g. "--max-size 200G"')
    ap.add_argument("--fasterq-extra", default="", help="Extra args for fasterq-dump")

    # metadata options
    ap.add_argument("--write-metadata", action="store_true", help="Write series matrix, runinfo, and merged metadata TSV")
    ap.add_argument("--metadata-only", action="store_true", help="Only write metadata and exit")

    # SLURM mode
    ap.add_argument("--slurm", action="store_true", help="Submit one SLURM job per run")
    ap.add_argument("--slurm-dry-run", action="store_true", help="Write sbatch scripts but do not submit")
    ap.add_argument("--slurm-jobs-dir", default=None, help="Where to write sbatch scripts (default: <outdir>/slurm_jobs)")
    ap.add_argument("--slurm-logdir", default=None, help="SLURM stdout/stderr dir (default: <outdir>/slurm_logs)")

    # Defaults requested
    ap.add_argument("--slurm-partition", default="standard", help="Default: standard")
    ap.add_argument("--slurm-account", default="dremel_lab", help="Default: dremel_lab")
    ap.add_argument("--slurm-time", default="12:00:00", help="Default: 12:00:00")
    ap.add_argument("--slurm-mem", default="48G", help="Default: 48G")

    ap.add_argument(
        "--slurm-prolog",
        default="source /project/dremel_lab/scripts/.sh_common\nmamba activate sra-tools",
        help="Shell snippet run in each SLURM job before commands. Default activates mamba env 'sra-tools'.",
    )

    args = ap.parse_args()

    gse = args.gse.strip()
    if not re.fullmatch(r"GSE\d+", gse):
        raise SystemExit(f"[error] Invalid GSE accession: {gse}")

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    pigz_threads = int(args.pigz_threads) if args.pigz_threads else int(args.threads)

    # metadata-only implies write-metadata
    if args.metadata_only:
        args.write_metadata = True

    if args.write_metadata:
        print(f"[info] Writing metadata for {gse} into {outdir}")
        write_metadata_bundle(gse, outdir)

    if args.metadata_only:
        print("[done] Metadata-only mode complete.")
        return

    runs = gse_to_runs(gse)
    if args.limit and args.limit > 0:
        runs = runs[: args.limit]

    prefetch_extra = args.prefetch_extra.split() if args.prefetch_extra.strip() else []
    fasterq_extra = args.fasterq_extra.split() if args.fasterq_extra.strip() else []

    if args.slurm:
        check_tool("sbatch")

        jobs_dir = os.path.abspath(args.slurm_jobs_dir or os.path.join(outdir, "slurm_jobs"))
        logdir = os.path.abspath(args.slurm_logdir or os.path.join(outdir, "slurm_logs"))

        print(f"[info] SLURM mode ON. Submitting {len(runs)} job(s).")
        if args.slurm_dry_run:
            print("[info] slurm dry-run: scripts will be written, not submitted.")

        for i, run_acc in enumerate(runs, 1):
            sbatch_script = write_sbatch_script(
                run_acc=run_acc,
                outdir=outdir,
                logdir=logdir,
                cpus=int(args.threads),
                mem=args.slurm_mem,
                time_limit=args.slurm_time,
                partition=args.slurm_partition,
                account=args.slurm_account,
                keep_sra=bool(args.keep_sra),
                pigz_on=bool(args.pigz),
                pigz_threads=pigz_threads,
                prefetch_extra=prefetch_extra,
                fasterq_extra=fasterq_extra,
                prolog=args.slurm_prolog,
            )
            submit_srr_as_slurm_job(
                run_acc=run_acc,
                jobs_dir=jobs_dir,
                sbatch_script=sbatch_script,
                dry_run=bool(args.slurm_dry_run),
            )
            print(f"[info] queued {i}/{len(runs)}: {run_acc}")

        print("[done] SLURM submissions complete.")
        return

    # Local mode
    check_tool("prefetch")
    check_tool("fasterq-dump")
    if args.pigz:
        check_tool("pigz")

    print(f"[info] Local mode. Processing {len(runs)} run(s) into {outdir}")
    for i, run_acc in enumerate(runs, 1):
        print(f"\n[{i}/{len(runs)}] {run_acc}")
        process_run_local(
            run_acc=run_acc,
            outdir=outdir,
            threads=int(args.threads),
            keep_sra=bool(args.keep_sra),
            pigz_on=bool(args.pigz),
            pigz_threads=pigz_threads,
            prefetch_extra=prefetch_extra,
            fasterq_extra=fasterq_extra,
        )
        if args.sleep:
            time.sleep(float(args.sleep))

    # Final sweep: remove any empty run directories left behind
    for name in os.listdir(outdir):
        if re.fullmatch(r"[SED]RR\d+", name):
            cleanup_empty_run_dir(outdir, name)

    print("\n[done] All runs processed (local mode).")


if __name__ == "__main__":
    main()

