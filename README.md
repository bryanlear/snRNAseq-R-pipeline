# snRNA-seq pipeline

This DAG shows the pipeline for the Alzheimer disease (AD) and control (Morabito et al., 2021 - [Nature Genetics](https://www.nature.com/articles/s41588-021-00894-z). Article can be found [here](literature/morabito_et_al_2021.pdf) and raw data [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367)).

1 `.rds` = 1 sample

```mermaid
---
config:
  themeVariables:
    fontSize: 12px
  flowchart:
    nodeSpacing: 20
    rankSpacing: 20
    padding: 8
    diagramPadding: 8
---

flowchart TB
    INPUT["AD and control data<br/>H5 count matrix + CSV metadata"]
    IMPORT["00_combined_h5.R<br/>Match barcodes and split by sample<br/>unmatched_policy = exclude"]
    IMPORT_OUT["Per-sample .rds files<br/>import_summary.tsv + excluded_barcodes.tsv"]
    QC["02_qc.R<br/>UMIs (nCount_RNA) and genes (nFeature_RNA)<br/>Filter log10 counts: median ± 3 × MAD"]
    QC_OUT["Filtered .rds + thresholds.tsv<br/>Count limits and removal summary"]
    MITO["03_mit_contamination.R<br/>Remove nuclei above mitochondrial cutoff<br/>Genes: ^MT-; method = mad<br/>Cutoff = median + max(3 × MAD, 0.5)"]
    MITO_OUT["Filtered .rds + summary.tsv"]
    CYCLE["04_cell_cycle_scoring.R<br/>LogNormalize scoring copy; scale.factor = 10000"]
    CYCLE_OUT[".rds + summary.tsv<br/>Append S.Score, G2M.Score, and Phase<br/>No regression"]
    DOUBLETS["05_doublets.R<br/>DoubletFinder; retain singlets<br/>PCs = 1:15; pN = 0.25<br/>Select pK at maximum BCmetric"]
    DOUBLETS_OUT["Singlet .rds + expected_doublets.tsv<br/>pK_selected.tsv + pK_sweep.tsv<br/>pK_sweep.png"]
    AMBIENT["06_ambient_rna.R<br/>decontX; seed = 12345<br/>z = NULL; background = NULL"]
    FINAL[".rds files<br/>Corrected decontX assay + contamination estimates"]

    INPUT --> IMPORT --> IMPORT_OUT
    IMPORT_OUT -->|.rds| QC --> QC_OUT
    QC_OUT -->|.rds| MITO --> MITO_OUT
    MITO_OUT -->|.rds| CYCLE --> CYCLE_OUT
    CYCLE_OUT -->|.rds| DOUBLETS --> DOUBLETS_OUT
    DOUBLETS_OUT -->|.rds| AMBIENT --> FINAL

    subgraph REPORTS["Optional QC reports — no nuclei removed"]
        QC_STATS["02.1_qc.R<br/>Mitochondrial gene pattern: ^MT-"]
        STATS_OUT["Full QC stats: cells.tsv + summary.tsv<br/>UMIs, genes, and mitochondrial RNA (%)"]
        QC_PLOTS["02.2_umi_gene_mito_plot.R<br/>nmads = 3; min_diff = 0.5"]
        PLOTS_OUT["3 .png plots + flagged.tsv<br/>Flag nuclei above the mitochondrial cutoff<br/>for removal review"]
        QC_STATS --> STATS_OUT
        STATS_OUT -->|cells.tsv| QC_PLOTS --> PLOTS_OUT
    end

    QC_OUT -.->|.rds| QC_STATS
    H5["Original H5 matrix<br/>Estimate recovered nuclei per sample"] --> DOUBLETS

    classDef process fill:#eaf2ff,stroke:#315a99,color:#172b4d;
    classDef data fill:#eaf7ee,stroke:#34734a,color:#173d25;
    classDef report fill:#fff4df,stroke:#a46b16,color:#583c13;
    class IMPORT,QC,MITO,CYCLE,DOUBLETS,AMBIENT process;
    class INPUT,IMPORT_OUT,QC_OUT,MITO_OUT,CYCLE_OUT,DOUBLETS_OUT,FINAL,H5 data;
    class QC_STATS,STATS_OUT,QC_PLOTS,PLOTS_OUT report;
```
