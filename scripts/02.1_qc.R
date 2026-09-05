args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L) {stop("Usage: Rscript scripts/02.1_qc.R <input.rds> <cells.tsv> <summary.tsv> <gene_pattern>", call. = FALSE)}

input_rds <- args[[1]]
cells_tsv <- args[[2]]
summary_tsv <- args[[3]]
mito_pattern <- args[[4]]

if (!file.exists(input_rds)) {stop("Input file doesn't exist: ", input_rds, call. = FALSE)}
if (!nzchar(mito_pattern)) {stop("Gene pattern can't be empty.", call. = FALSE)}

suppressPackageStartupMessages(library(Seurat))

object <- readRDS(input_rds)

if (!inherits(object, "Seurat") || ncol(object) == 0L) {stop("Input must be a Seurat object with nuclei.", call. = FALSE)}
if (!"RNA" %in% Assays(object)) {stop("RNA assay is required.", call. = FALSE)}

metadata <- object[[]]
required <- c("sample_id", "condition", "nCount_RNA", "nFeature_RNA")

if (!all(required %in% names(metadata))) {stop("Metadata requires sample_id, condition, nCount_RNA, and nFeature_RNA.", call. = FALSE)}

sample_id <- unique(as.character(metadata$sample_id))
condition <- unique(as.character(metadata$condition))

if (length(sample_id) != 1L || anyNA(sample_id) || !nzchar(sample_id)) {stop("Input must contain 1 nonempty sample_id.", call. = FALSE)}
if (length(condition) != 1L || anyNA(condition) || !nzchar(condition)) {stop("Input must contain 1 nonempty condition.", call. = FALSE)}

counts_layers <- Layers(object[["RNA"]], search = "^counts($|\\.)")

if (length(counts_layers) != 1L) {stop("RNA assay must contain exactly one counts layer. Join split counts layers before this step.", call. = FALSE)}

counts <- GetAssayData(object, assay = "RNA", layer = counts_layers[[1]])

if (!identical(colnames(counts), colnames(object)) || !identical(rownames(metadata), colnames(object))) {stop("RNA counts and metadata must contain the same nuclei in the same order.", call. = FALSE)}

total_umis <- Matrix::colSums(counts)
detected_genes <- Matrix::colSums(counts > 0)

if (any(!is.finite(total_umis)) || any(total_umis <= 0)) {stop("Each nucleus must have a positive, finite total UMI count.", call. = FALSE)}
if (!isTRUE(all.equal(unname(total_umis), as.numeric(metadata$nCount_RNA))) || !isTRUE(all.equal(unname(detected_genes), as.numeric(metadata$nFeature_RNA)))) {stop("RNA counts do not agree with nCount_RNA or nFeature_RNA metadata.", call. = FALSE)}

mito_genes <- grep(mito_pattern, rownames(counts), value = TRUE)

if (length(mito_genes) == 0L) {stop("No mitochondrial genes match the pattern.", call. = FALSE)}

mito_counts <- counts[mito_genes, , drop = FALSE]
mitochondrial_umis <- Matrix::colSums(mito_counts)
detected_mitochondrial_genes <- Matrix::colSums(mito_counts > 0)
percent_mt <- 100 * mitochondrial_umis / total_umis

if (anyNA(percent_mt) || any(!is.finite(percent_mt)) || any(percent_mt < 0 | percent_mt > 100)) {stop("Mitochondrial percentages must be finite values between 0 and 100.", call. = FALSE)}

report_columns <- c(required, intersect(c("Batch", "Cell.Type", "cluster"), names(metadata)))
cell_report <- data.frame(Barcode = colnames(object), metadata[, report_columns, drop = FALSE], mitochondrial_umis = unname(mitochondrial_umis), detected_mitochondrial_genes = unname(detected_mitochondrial_genes), percent.mt = unname(percent_mt), check.names = FALSE)

q <- quantile(percent_mt, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1), names = FALSE)

summary <- data.frame(
  sample_id = sample_id, condition = condition, nuclei = ncol(object), mito_pattern = mito_pattern, mitochondrial_features = length(mito_genes),
  total_umis = sum(total_umis), total_mitochondrial_umis = sum(mitochondrial_umis), genes_detected_in_sample = sum(Matrix::rowSums(counts > 0) > 0), mitochondrial_genes_detected_in_sample = sum(Matrix::rowSums(mito_counts > 0) > 0),
  min_umis_per_nucleus = min(total_umis), median_umis_per_nucleus = median(total_umis), mean_umis_per_nucleus = mean(total_umis), max_umis_per_nucleus = max(total_umis),
  min_genes_per_nucleus = min(detected_genes), median_genes_per_nucleus = median(detected_genes), mean_genes_per_nucleus = mean(detected_genes), max_genes_per_nucleus = max(detected_genes),
  min_mt_umis_per_nucleus = min(mitochondrial_umis), median_mt_umis_per_nucleus = median(mitochondrial_umis), mean_mt_umis_per_nucleus = mean(mitochondrial_umis), max_mt_umis_per_nucleus = max(mitochondrial_umis),
  pooled_percent_mt = 100 * sum(mitochondrial_umis) / sum(total_umis),
  nuclei_zero_mt = sum(percent_mt == 0), percent_nuclei_zero_mt = 100 * mean(percent_mt == 0),
  min_percent_mt = q[1], q25_percent_mt = q[2], median_percent_mt = q[3], mean_percent_mt = mean(percent_mt), mad_percent_mt = mad(percent_mt),
  q75_percent_mt = q[4], q90_percent_mt = q[5], q95_percent_mt = q[6], q99_percent_mt = q[7], max_percent_mt = q[8])

dir.create(dirname(cells_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_tsv), recursive = TRUE, showWarnings = FALSE)

write.table(cell_report, file = cells_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file = summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

message("Checked ", sample_id, ": ", ncol(object), " nuclei; ", sum(total_umis), " total UMIs; ", sum(mitochondrial_umis), " mitochondrial UMIs; ", summary$genes_detected_in_sample, " detected genes; median mitochondrial RNA = ", signif(median(percent_mt), 6), "%; maximum = ", signif(max(percent_mt), 6), "%.")
