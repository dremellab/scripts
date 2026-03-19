#!/usr/bin/env python3
"""
transfer_pipeline_results.py

Create Globus batch files for pipeline outputs based on lookup.tsv.
Each pipeline has its own transfer rules (extensions + destination layout).

Examples:
  # dry run to show transfer command
  transfer_pipeline_results.py --sampleSetName 20250423_RNAseq-1 --dry-run

  # build filelist/batchfile and submit transfer
  transfer_pipeline_results.py --sampleSetName 20250423_RNAseq-1
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
from typing import Dict, Iterable, List, Tuple


DEFAULT_PROJECT = os.environ.get("PROJECT", "/project/dremel_lab")
DEFAULT_LOOKUP = os.path.join(DEFAULT_PROJECT, "analysis", "lookup.tsv")
DEFAULT_ANALYSIS_DIR = os.path.join(DEFAULT_PROJECT, "analysis")
DEFAULT_SOURCE_ROOT = "/dtn/landings/users/c/cu/cud2td/project/dremel_lab/analysis"
DEFAULT_SOURCE_UUID = "af187d15-768f-4449-8670-d00e1eb1ce6a"
DEFAULT_DESTINATION_UUID = "af187d15-768f-4449-8670-d00e1eb1ce6a"
LAPTOP_DESTINATION_UUID = "6bfd96d9-050c-11f0-ad00-0e283342ad7b"
LAPTOP_DEST_ROOT = "/Users/vishal/Documents/Data/Analysis/{sampleSetName}"
S3_DESTINATION_UUID = "577d6907-4263-49a4-b0c7-3f6b80064d0b"
S3_DEST_ROOT = "/dremel-lab-bucket/_HTS/{sampleSetName}"


BAM_RULE = {
    "kind": "suffix",
    "suffixes": [".bam", ".bai"],
    "dest": "bams",
    "exclude_substrings": ["Aligned.out.bam"],
}

PIPELINE_CONFIG: Dict[str, Dict] = {
    "harold": {
        "dest_root": "/dtn/landings/storage/leased/vol_dremellab/_HTS/{sampleSetName}/_Outputs",
        "rules": [
            {"kind": "path", "path": "samples.tsv", "dest": "config/samples.tsv"},
            {"kind": "path", "path": "config.yaml", "dest": "config/config.yaml"},
            {"kind": "path", "path": "config/rivanna/config.yaml", "dest": "config/rivanna/config.yaml"},
            {"kind": "path", "path": "results/alignmentqc/alignment_summary.tsv", "dest": "results/alignmentqc/alignment_summary.tsv"},
            {"kind": "suffix", "suffixes": [".bw", ".bb"], "dest": "bigwigs"},
            {"kind": "dir", "dir_name": "counts", "dest": "counts"},
            {"kind": "dir", "dir_name": "multiqc_data", "dest": "multiqc_data"},
            {"kind": "basename", "basename": "multiqc_report.html", "dest": ""},
        ],
    },
}


def load_lookup_table(path: str) -> Dict[str, Dict[str, str]]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Lookup table not found: {path}")
    rows: Dict[str, Dict[str, str]] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "sampleSetName" not in reader.fieldnames or "pipelineName" not in reader.fieldnames:
            raise ValueError("lookup.tsv must include sampleSetName and pipelineName columns")
        for row in reader:
            sample = (row.get("sampleSetName") or "").strip()
            if not sample:
                continue
            if sample in rows:
                raise ValueError(f"Duplicate sampleSetName in lookup.tsv: {sample}")
            rows[sample] = row
    return rows


def resolve_workdir(row: Dict[str, str], sample_set: str, analysis_dir: str) -> str:
    workdir = (row.get("workdir") or "").strip()
    if workdir:
        if os.path.isabs(workdir):
            return workdir
        return os.path.join(analysis_dir, workdir)
    return os.path.join(analysis_dir, sample_set)


def gather_files(source_dir: str) -> List[str]:
    relpaths: List[str] = []
    for root, dirs, files in os.walk(source_dir):
        dirs.sort()
        files.sort()
        for name in files:
            full = os.path.join(root, name)
            rel = os.path.relpath(full, source_dir)
            relpaths.append(rel.replace(os.sep, "/"))
    return relpaths


def write_filelist(path: str, relpaths: Iterable[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for rel in relpaths:
            f.write(f"{rel}\n")


def read_filelist(path: str) -> List[str]:
    relpaths: List[str] = []
    with open(path) as f:
        for line in f:
            rel = line.strip()
            if not rel or rel.endswith("/"):
                continue
            relpaths.append(rel)
    return relpaths


def match_rule(relpath: str, rule: Dict) -> str | None:
    kind = rule["kind"]
    if kind == "path":
        if relpath == rule.get("path"):
            return rule.get("dest")
    elif kind == "suffix":
        for suffix in rule.get("suffixes", []):
            if relpath.endswith(suffix):
                if any(s in relpath for s in rule.get("exclude_substrings", [])):
                    return None
                dest_dir = rule.get("dest", "")
                base = os.path.basename(relpath)
                return f"{dest_dir}/{base}" if dest_dir else base
    elif kind == "basename":
        if os.path.basename(relpath) == rule.get("basename"):
            dest_dir = rule.get("dest", "")
            base = os.path.basename(relpath)
            return f"{dest_dir}/{base}" if dest_dir else base
    elif kind == "dir":
        dir_name = rule.get("dir_name")
        if not dir_name:
            return None
        parts = relpath.split("/")
        if dir_name in parts:
            idx = parts.index(dir_name)
            subpath = "/".join(parts[idx + 1 :])
            if not subpath:
                return None
            dest_dir = rule.get("dest", dir_name)
            return f"{dest_dir}/{subpath}"
    return None


def build_batch_entries(relpaths: Iterable[str], rules: List[Dict]) -> List[Tuple[str, str]]:
    entries: List[Tuple[str, str]] = []
    seen: Dict[str, str] = {}
    for rel in relpaths:
        for rule in rules:
            dest = match_rule(rel, rule)
            if dest:
                if rel in seen and seen[rel] != dest:
                    raise ValueError(f"Conflicting destinations for {rel}: {seen[rel]} vs {dest}")
                if rel not in seen:
                    seen[rel] = dest
                    entries.append((rel, dest))
                break
    return entries


def ensure_trailing_slash(path: str) -> str:
    return path if path.endswith("/") else path + "/"


def ensure_globus_env() -> None:
    if os.environ.get("TRANSFER_PIPELINE_RESULTS_ENV_CHECK") == "1":
        return
    current_env = os.environ.get("CONDA_DEFAULT_ENV") or os.environ.get("MAMBA_DEFAULT_ENV")
    if current_env == "globus":
        return
    script = os.path.abspath(sys.argv[0])
    args = " ".join(shlex.quote(arg) for arg in sys.argv[1:])
    cmd = (
        "source /project/dremel_lab/scripts/.sh_common && "
        "mamba activate globus && "
        f"TRANSFER_PIPELINE_RESULTS_ENV_CHECK=1 {shlex.quote(sys.executable)} {shlex.quote(script)} {args}"
    )
    result = subprocess.run(cmd, shell=True, executable="/bin/bash")
    if result.returncode != 0:
        raise SystemExit("Failed to activate mamba env 'globus'.")
    raise SystemExit(result.returncode)


def main() -> None:
    ensure_globus_env()
    parser = argparse.ArgumentParser(description="Create Globus batch files for pipeline results.")
    parser.add_argument("--sampleSetName", required=True, help="Sample set name (lookup.tsv key)")
    parser.add_argument("--lookup", default=DEFAULT_LOOKUP, help="Path to lookup.tsv")
    parser.add_argument("--analysis-dir", default=DEFAULT_ANALYSIS_DIR, help="Local analysis root")
    parser.add_argument("--source-root", default=DEFAULT_SOURCE_ROOT, help="Globus source root")
    parser.add_argument(
        "--transfer-type",
        choices=("rivanna-to-topaz", "rivanna-to-laptop", "rivanna-to-s3"),
        default="rivanna-to-topaz",
        help="Transfer route to select default endpoint UUIDs",
    )
    parser.add_argument(
        "--source-uuid",
        default=DEFAULT_SOURCE_UUID,
        help="Globus source endpoint UUID",
    )
    parser.add_argument(
        "--destination-uuid",
        default=DEFAULT_DESTINATION_UUID,
        help="Globus destination endpoint UUID",
    )
    parser.add_argument("--filelist", help="Override filelist path")
    parser.add_argument("--batchfile", help="Override batchfile path")
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="Keep generated filelist and batchfile (default: delete after transfer)",
    )
    parser.add_argument(
        "--include-bam",
        action="store_true",
        help="Include .bam and .bai files in transfer (excluded by default)",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print transfer command only")
    args = parser.parse_args()
    if args.transfer_type == "rivanna-to-laptop" and args.destination_uuid == DEFAULT_DESTINATION_UUID:
        args.destination_uuid = LAPTOP_DESTINATION_UUID
    elif args.transfer_type == "rivanna-to-s3" and args.destination_uuid == DEFAULT_DESTINATION_UUID:
        args.destination_uuid = S3_DESTINATION_UUID

    lookup = load_lookup_table(args.lookup)
    row = lookup.get(args.sampleSetName)
    if not row:
        raise SystemExit(f"SampleSetName not found in lookup.tsv: {args.sampleSetName}")
    pipeline = (row.get("pipelineName") or "").strip()
    if not pipeline:
        raise SystemExit(f"pipelineName missing for sampleSetName: {args.sampleSetName}")

    config = PIPELINE_CONFIG.get(pipeline)
    if not config:
        raise SystemExit(f"No transfer config for pipeline: {pipeline}")
    rules = config.get("rules", [])
    if not rules:
        raise SystemExit(f"No transfer rules configured for pipeline: {pipeline}")
    
    # Conditionally add BAM rule if --include-bam is specified
    if args.include_bam:
        rules = rules + [BAM_RULE]

    analysis_dir = args.analysis_dir
    source_dir = resolve_workdir(row, args.sampleSetName, analysis_dir)
    if not os.path.isdir(source_dir):
        raise SystemExit(f"Source directory not found: {source_dir}")

    filelist_path = args.filelist or os.path.join(analysis_dir, f"{args.sampleSetName}.filelist")
    batchfile_path = args.batchfile or os.path.join(analysis_dir, f"{args.sampleSetName}.batchfile")

    relpaths = gather_files(source_dir)
    write_filelist(filelist_path, relpaths)
    print(f"[filelist] wrote {len(relpaths)} entries to {filelist_path}")

    relpaths = read_filelist(filelist_path)
    entries = build_batch_entries(relpaths, rules)
    if not entries:
        raise SystemExit("No files matched transfer rules.")
    os.makedirs(os.path.dirname(batchfile_path), exist_ok=True)
    with open(batchfile_path, "w") as f:
        for src, dst in entries:
            f.write(f"{src} {dst}\n")
    print(f"[batchfile] wrote {len(entries)} entries to {batchfile_path}")

    workdir = resolve_workdir(row, args.sampleSetName, analysis_dir)
    if os.path.isabs(workdir):
        try:
            rel_workdir = os.path.relpath(workdir, analysis_dir)
        except ValueError:
            rel_workdir = None
    else:
        rel_workdir = workdir
    if not rel_workdir or rel_workdir.startswith(".."):
        raise SystemExit(
            "Cannot map absolute workdir outside analysis-dir to globus source-root. "
            "Provide a compatible --analysis-dir/--source-root or update lookup.tsv."
        )

    source_path = ensure_trailing_slash(os.path.join(args.source_root, rel_workdir))
    if args.transfer_type == "rivanna-to-laptop":
        dest_root = LAPTOP_DEST_ROOT.format(
            sampleSetName=args.sampleSetName,
            pipelineName=pipeline,
        )
    elif args.transfer_type == "rivanna-to-s3":
        dest_root = S3_DEST_ROOT.format(
            sampleSetName=args.sampleSetName,
            pipelineName=pipeline,
        )
    else:
        dest_root = config["dest_root"].format(
            sampleSetName=args.sampleSetName,
            pipelineName=pipeline,
        )
    dest_path = ensure_trailing_slash(dest_root)
    cmd = [
        "globus",
        "transfer",
        "--batch",
        batchfile_path,
        f"{args.source_uuid}:{source_path}",
        f"{args.destination_uuid}:{dest_path}",
    ]
    if args.dry_run:
        print("[dry-run]", " ".join(cmd))
        return
    subprocess.run(cmd, check=True)
    if not args.keep_intermediate:
        for path in (filelist_path, batchfile_path):
            try:
                os.remove(path)
            except FileNotFoundError:
                pass


if __name__ == "__main__":
    main()
