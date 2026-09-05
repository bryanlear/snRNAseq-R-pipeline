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

mito_genes <- grep(mito_pattern, rownames(object[["RNA"]]), value = TRUE)

if (length(mito_genes) == 0L) {stop("No mitochondrial genes match the pattern.", call. = FALSE)}

percent_mt <- PercentageFeatureSet(object, features = mito_genes, assay = "RNA")

if (anyNA(percent_mt) || any(!is.finite(percent_mt)) || any(percent_mt < 0 | percent_mt > 100)) {stop("Mitochondrial percentages must be finite values between 0 and 100.", call. = FALSE)}

report_columns <- c(required, intersect(c("Batch", "Cell.Type", "cluster"), names(metadata)))
cell_report <- data.frame(Barcode = colnames(object), metadata[, report_columns, drop = FALSE], percent.mt = unname(percent_mt), check.names = FALSE)

q <- quantile(percent_mt, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1), names = FALSE)

summary <- data.frame(
  sample_id = sample_id, condition = condition, nuclei = ncol(object), mito_pattern = mito_pattern, mitochondrial_features = length(mito_genes),
  nuclei_zero_mt = sum(percent_mt == 0), percent_nuclei_zero_mt = 100 * mean(percent_mt == 0),
  min_percent_mt = q[1], q25_percent_mt = q[2], median_percent_mt = q[3], mean_percent_mt = mean(percent_mt), mad_percent_mt = mad(percent_mt),
  q75_percent_mt = q[4], q90_percent_mt = q[5], q95_percent_mt = q[6], q99_percent_mt = q[7], max_percent_mt = q[8])

dir.create(dirname(cells_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_tsv), recursive = TRUE, showWarnings = FALSE)

write.table(cell_report, file = cells_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file = summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

message("Checked ", sample_id, ": ", ncol(object), " nuclei; ", length(mito_genes), " mitochondrial features; median = ", signif(median(percent_mt), 6), "%; maximum = ", signif(max(percent_mt), 6), "%.")
