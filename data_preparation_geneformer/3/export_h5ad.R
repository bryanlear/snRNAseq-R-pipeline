### Export `sample_id`, `condition`, `batch`, `cell.type` data from geneformer_rds dir to .h5ad for Geneformer input

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

input_dir <- "data_preparation_geneformer/geneformer_rds"
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)

stopifnot(length(rds_files) > 0L)

for (input_file in rds_files) {object <- readRDS(input_file)
    for (assay_name in c("RNA", "decontX")){counts <- LayerData(object, assay = assay_name, layer = "counts")

        cell_metadata <- object[[]][colnames(counts), c("sample_id", "condition", "Batch", "Cell.Type"), drop = FALSE]
        colnames(cell_metadata) <- c("sample_id", "condition", "batch", "cell.type")

        gene_metadata <- object[[assay_name]][[]][rownames(counts), , drop = FALSE]

        sce <- SingleCellExperiment(assays = list(counts = counts), colData = S4Vectors::DataFrame(cell_metadata), rowData = S4Vectors::DataFrame(gene_metadata))

        output_dir <- file.path("data_preparation_geneformer", paste0("h5ad_", assay_name))
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        output_file <- file.path(output_dir, sub("\\.rds$", ".h5ad", basename(input_file)))

        writeH5AD(sce, file = output_file, X_name = "counts", compression = "gzip")
        message("Shit saved: ", output_file)}

    rm(object, counts, sce)
    gc()}