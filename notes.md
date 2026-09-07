
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

6. **Append Cell-Cycle Scoring Medata to Object**

```
snakemake cell_cycle_scoring_all --dry-run --printshellcmds --cores 1
```

7. **Doublets**:

First PCA and investigate 8 to 15 PCs:

e.g., Sample 22
```
Rscript scripts/05_doublets.R \
    results/GSE174367/04_cell_cycle/Sample-22.rds \
    results/GSE174367/05_doublets/Sample-22.rds \
    results/GSE174367/05_doublets/Sample-22.summary.tsv

---

PC_ 8 
Positive:  RALYL, SH3GL2, UNC13C, GAP43, AC073365.1 
Negative:  GALNTL6, CALB2, TAC3, BTBD11, PWRN1 
PC_ 9 
Positive:  PTPRM, NEAT1, MYRIP, AFF3, RASGRF1 
Negative:  ADGRL3, AL445250.1, CTNND2, KCNAB1, SVEP1 
PC_ 10 
Positive:  TLL1, SPHKAP, RIT2, EBF1, ZNF385D 
Negative:  SYT7, GYG2P1, TMEM196, MLIP, AC073365.1 
PC_ 11 
Positive:  RAB3C, SATB1-AS1, TSHZ3, ST6GALNAC5, BCL11B 
Negative:  PPFIA4, LRRC3B, EGFEM1P, GPR158, LINGO2 
PC_ 12 
Positive:  CCDC26, IRF8, MINOS1, KCNIP1, OR52K3P 
Negative:  CD163, IQGAP2, F13A1, SLC16A10, MRC1 
PC_ 13 
Positive:  SGCZ, PTPRT, HIF1A-AS2, SYTL3, GPNMB 
Negative:  VSIG4, STAB1, CD163, MRC1, ANKRD55 
PC_ 14 
Positive:  ANKRD55, EGFEM1P, SPTLC3, DACH2, BEST3 
Negative:  SGCZ, VSIG4, GPNMB, NCK2, STAB1 
PC_ 15 
Positive:  LYPD6, RERG, NPNT, CHI3L1, PPFIA4 
Negative:  HMGCLL1, ABI3BP, ANGPT1, FRMPD2, PRR16 
```
[Source DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder#doubletfinder-overview)

Note: Positive and negative loadings indicate direction along PC and NOT differential expression. 

* **pK**: Defines PC neighborhood size/radius used to compute proportion of artificial nearest neighbors (**pANN**) for each cell. It is the radius as a fracation of the total merged real and artificial dataset. Since optimal performance depends of dataset density and cell-type composition, pK, is thus not assigned a default value. Use `paramSweep` and `find.pK` to find peak in mean variance normalized bimodality coefficient (**BCmvn**).