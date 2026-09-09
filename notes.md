# PRE-PROCESSING:

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
    results/GSE174367/05_doublets/Sample-22.summary.tsv \
    brain_data_human/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5 --> Use these counts to calculate n_recovered and expected rate
```

[Source DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder#doubletfinder-overview)

Note: Positive and negative loadings indicate direction along PC and NOT differential expression. 

* **pK**: Defines PC neighborhood size/radius used to compute proportion of artificial nearest neighbors (**pANN**) for each cell. It is the radius as a fracation of the total merged real and artificial dataset. Since optimal performance depends of dataset density and cell-type composition, pK, is thus not assigned a default value. Use `paramSweep` and `find.pK` to find peak in mean variance normalized bimodality coefficient (**BCmvn**).

```
snakemake doublet_finder_all --cores 3 --printshellcmds --allowed-rules doublet_finder doublet_finder_all
```

8. **Ambient RNA**:

**DecontX**: Bayesian hirarchical method which estimates and remove ambient RNA contaimination (a.k.a *soup*). Ambient transcripts which are released from lyzed cells during tissue dissociation, are captured inside the droplets along the cells/nuclei native mRNA. DecontX decouples those signals withouth requiring empty-droplet information. 

DecontX assumes that ambient RNA pool (soup) consists of a uniform mizture of transcripts leaked from all lyzed cells in the suspension. 

$$P(X_{ij} \mid \psi_i, \phi_{k[i], j}, \gamma_j) = (1 - \psi_i)\phi_{k[i], j} + \psi_i \gamma_j$$

 where $\gamma_j$ is global background freq. of gene $j$ in the ambient RNA pool.

```
snakemake ambient_rna_all --cores 3 --printshellcmds --allowed-rules ambient_rna ambient_rna_all
```

**NOTE regarding ambient RNA**: I checked correction using UMI changes and whatnot. They show the amount of correction but do NOT establish biological accuracy. **For a future research project**: Compare several cell-type marker genes before and after correction (final saved `.rds` since object contains both columns). Must check that expected expression remains in the appropriate cell types and unexpected decreases somewhere else. Must review nuclei/cell with large corrections and check cell groups used by the tool (here `DecontX`). Any filtering step must have justification. 

**UMAP**:





---

# GENEFORMER v2 (104 * 10^6 Human Transcriptomes)

Chen et al,. 2026 [Nature Computational Science](https://www.nature.com/articles/s43588-026-00972-4). Article can be found [here](literature/chen_et_al_2026.pdf).

Counts for each nucleus $\rightarrow$ Gene IDs and metadata $\rightarrow$ Tokenizer: calculate scores and rank genes $\rightarrow$ Geneformer: process ranked gene tokens $\rightarrow$ Embeddings $\rightarrow$ Analysis

Tokenizer calculates a score for each detected gene:

$$s_{g,c}=\frac{10,000x_{g,c}}{N_c m_g}$$

- $x_{g,c}$ is count for gene $g$ in nucleus $c$
- $N_c$ is total count in $c$
- $m_g$ is reference median for $g$

Tokenizer sorts genes from highest score to lowest. Reference medians come from pretraining corpus.

## Data Preparation:

1. Use DecontX counts $\rightarrow$ compare corrected vs. non-corrected which are both in the .rds
2. Restore Ensembl IDs (from original .H5)
3. Export to `.h5ad` $\rightarrow$ `sample_id`, `condition`, `batch`, `cell.type`
4. Tokenize $\rightarrow$ feed counts with NO log tansform, scaling, and whatnot. Tokenizer takes care of that.

### 1.

**Pericyte cells**: Mural cells embedded within the vascular basement membrane that wrap around endothelial cells lining capillaries and post capillary venules. They have roles in vascular stability/integrity, BBB regulation, microvascular blood flow, angiogenesis.

`PER.END` = Pericytes and Endothelial cells annotation in the `.rsd` files.
