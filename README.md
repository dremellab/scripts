## Dremel Lab common scripts folder

- Location on Rivanna: `/project/dremel_lab/scripts`
- All `dremel_lab` group members should add the following to their `~/.bashrc` files. In addition to making all shared scripts like `jobby` available in path, this also makes shared `mamba` and its environments available to all. For example, now you can run `mamba activate pipelines` to have all pipelines ready to go!! 

```bash
source /project/dremel_lab/scripts/.sh_common
```

