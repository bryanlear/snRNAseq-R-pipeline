`commandArgs()` = (with `trailingOnly = FALSE` as default) returns all arguments passed to the R process, including path to R executable itself and standard R invocation flags (e.g., `--no-echo`, `--vanilla`, `-f`):
    - `--no-echo`: Disables echoing of input commands to stdout
    - `--vanilla`: Forces R to start a clean, lightweight session with default setting

With `trailingOnly = TRUE`, it returns only the custom arguments supplied after the script name when running via `Rscript` or `R --slave -f`.

1. **TEST: Create Seurat Object:**

```
Rscript scripts/01_create_seurat.R \
  [sample_id] \
  [experimental_condition] \
  [file.h5] \
  results/01_seurat/output_files.rds
```