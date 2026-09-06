
1. **TEST: Handle Human Alzheimer's Data + Metada**

```
Rscript scripts/00_combined_h5.R \
  brain_data_human/snRNA-matrix.h5 \
  brain_data_human/meta.csv.gz \
  results/GSE174367/01_seurat \
  exclude
```

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

---

***QC For Human data:***

```
snakemake human_all --cores 1 --printshellcmds
````
---

3. **Mitochondrial contamination percentages summary**

```
snakemake human_mito_check_all --printshellcmds --cores 1
```

4. **Quality Control Plots**

```
snakemake human_qc_plots_all --cores 1 --printshellcmds
```

-  Also added to the plots an upper threshold = median_mt + 3 (multiplier) * MAD_mt and nuclei above is flagged. 
-  No nuclei is removed, only for visualization and inspection of the nuclei with relatively high mitochondrial RNA.

For context:
```median_mt <- median(percent_mt)

distances <- abs(percent_mt - median_mt)

mad_mt <- 1.4826 * median(distances)  # same as R's mad(percent_mt)

upper_threshold <- median_mt + 3 * mad_mt
```
5. **Apply automated mitochondrial decontamination thresholds**

```
snakemake mit_decontamination_all --printshellcmds --cores 1
```