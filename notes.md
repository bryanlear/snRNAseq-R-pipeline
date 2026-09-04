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

2. **TEST: QC (UMIs and Gene count using MADs)**

```
log_values <- log10(values)
center <- median(log_values)
spread <- mad(log_values)

and then:

lower <- 10^(center - 3 * spread)
upper <- 10^(center + 3 * spread)

and it is applied separately to:

nCount_RNA and nFeature_RNA
```

```
Rscript scripts/02_qc.R \
  results/01_seurat/output_files.rds \
  results/02_qc/output_files.rds \
  results/02_qc/output_files.thresholds.tsv \
  3 
```