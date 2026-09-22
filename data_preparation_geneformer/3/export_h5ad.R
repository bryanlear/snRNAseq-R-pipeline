### Export raw counts, Geneformer-required metadata, and cell metadata to .h5ad.

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

input_dir <- "data_preparation_geneformer/geneformer_rds"
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

stopifnot(length(rds_files) > 0L)

for (input_file in rds_files) {object <- readRDS(input_file)

    ensembl_metadata <- object[["decontX"]][[]]

    stopifnot(all(c("ensembl_id", "gene_symbol") %in% colnames(ensembl_metadata)), all(grepl("^ENSG[0-9]+$", ensembl_metadata$ensembl_id)))

    ensembl_ids <- as.character(ensembl_metadata$ensembl_id)
    rna_feature_names <- make.unique(as.character(ensembl_metadata$gene_symbol))

    for (assay_name in c("RNA", "decontX")) {counts <- LayerData(object, assay = assay_name, layer = "counts")

        stopifnot(nrow(counts) == length(ensembl_ids), all(counts@x >= 0), all(counts@x == round(counts@x)))

        if (assay_name == "RNA") {stopifnot(identical(rownames(counts), rna_feature_names))} else {stopifnot(identical(rownames(counts), rownames(ensembl_metadata)))}

        cell_metadata <- object[[]][colnames(counts), c("sample_id", "condition", "Batch", "Cell.Type"), drop = FALSE]
        colnames(cell_metadata) <- c("sample_id", "condition", "batch", "cell.type")
        cell_metadata$n_counts <- as.numeric(Matrix::colSums(counts))

        stopifnot(identical(rownames(cell_metadata), colnames(counts)), all(is.finite(cell_metadata$n_counts)), all(cell_metadata$n_counts > 0))

        gene_metadata <- object[[assay_name]][[]][rownames(counts), , drop = FALSE]
        gene_metadata$ensembl_id <- ensembl_ids

        sce <- SingleCellExperiment(assays = list(counts = counts), colData = S4Vectors::DataFrame(cell_metadata), rowData = S4Vectors::DataFrame(gene_metadata))

        output_dir <- file.path("data_preparation_geneformer", paste0("h5ad_", assay_name))
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        output_file <- file.path(output_dir, sub("\\.rds$", ".h5ad", basename(input_file)))

        writeH5AD(sce, file = output_file, X_name = "counts",compression = "gzip")
        message("Saved: ", output_file)}

    rm(object, counts, sce)
    gc()}
