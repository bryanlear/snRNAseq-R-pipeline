args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L){stop("Usage: Rscript scripts/04_cell_cycle_scoring.R <input.rds> <output.rds> <summary.tsv>", call. = FALSE)}

input_rds <- args[[1]]
output_rds <- args[[2]]
summary_tsv <- args[[3]]

if (!file.exists(input_rds)){stop("File doesn't exist: ", input_rds, call. = FALSE)}
if (anyDuplicated(normalizePath(c(input_rds, output_rds, summary_tsv), mustWork = FALSE))){stop("Use different paths for input and outputs.", call. = FALSE)}

suppressPackageStartupMessages(library(Seurat))
object <- readRDS(input_rds)

if (!inherits(object, "Seurat") || ncol(object) == 0L) {stop("Input must be Seurat object containing nuclei.", call. = FALSE)}
if (!"RNA" %in% Assays(object)) {stop("RNA assay is needed.", call. = FALSE)}

metadata <- object[[]]

if (!all(c("sample_id", "condition") %in% names(metadata))) {stop("Metadata requires sample_id and condition.", call. = FALSE)}
if (any(c("S.Score", "G2M.Score", "Phase") %in% names(metadata))) {stop("Cell-cycle metadata already exists. Use object from previous pipeline step.", call. = FALSE)}

sample_id <- unique(as.character(metadata$sample_id))
condition <- unique(as.character(metadata$condition))

if (length(sample_id) != 1L || anyNA(sample_id) || !nzchar(sample_id)) {stop("Input must contain 1 nonempty sample_id.", call. = FALSE)}
if (length(condition) != 1L || anyNA(condition) || !nzchar(condition)) {stop("Input must contain 1 nonempty condition.", call. = FALSE)}

counts_layers <- Layers(object[["RNA"]], search = "^counts($|\\.)")

if (length(counts_layers) != 1L || counts_layers[[1]] != "counts") {stop("RNA assay requires 1 counts layer named counts. Join split counts layers before this step.", call. = FALSE)}

counts <- GetAssayData(object, assay = "RNA", layer = "counts")

if (!identical(colnames(counts), colnames(object))) {stop("RNA counts must contain all nuclei in the same order as the object.", call. = FALSE)}
if (any(!is.finite(Matrix::colSums(counts))) || any(Matrix::colSums(counts) <= 0)) {stop("Each nucleus must have a positive, finite total count.", call. = FALSE)}

#for human only using 2019 Seurat's list
s_markers <- Seurat::cc.genes.updated.2019$s.genes
g2m_markers <- Seurat::cc.genes.updated.2019$g2m.genes
s_genes <- intersect(s_markers, rownames(counts))
g2m_genes <- intersect(g2m_markers, rownames(counts))
s_missing <- setdiff(s_markers, s_genes)
g2m_missing <- setdiff(g2m_markers, g2m_genes)

if (length(s_genes) == 0L || length(g2m_genes) == 0L) {stop("Both cell-cycle marker sets require matching human gene symbols.", call. = FALSE)}
if (length(s_missing) > 0L || length(g2m_missing) > 0L) {warning("Missing S markers: ", paste(s_missing, collapse = ","), "; missing G2M markers: ", paste(g2m_missing, collapse = ","), call. = FALSE)}

s_detected <- sum(Matrix::rowSums(counts[s_genes, , drop = FALSE]) > 0)
g2m_detected <- sum(Matrix::rowSums(counts[g2m_genes, , drop = FALSE]) > 0)

if (s_detected == 0L || g2m_detected == 0L) {stop("Both marker sets require detected expression for scoring.", call. = FALSE)}

# Normalize only scoring copy. NO SCALING/REGRESSION
scoring_object <- NormalizeData(object, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

#done on copy
scoring_object <- CellCycleScoring(scoring_object, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE, assay = "RNA", slot = "data", nbin = 24, seed = 1, search = FALSE)

#append to original
scores <- scoring_object[[]][colnames(object), c("S.Score", "G2M.Score", "Phase"), drop = FALSE]

if (any(!is.finite(scores$S.Score)) || any(!is.finite(scores$G2M.Score)) || anyNA(scores$Phase)) {stop("Cell-cycle scoring produced missing or invalid values.", call. = FALSE)}

object <- AddMetaData(object, metadata = scores)
rm(scoring_object)

phase_counts <- table(factor(scores$Phase, levels = c("G1", "S", "G2M", "Undecided")))

summary <- data.frame(
    sample_id = sample_id, condition = condition, nuclei = ncol(object), marker_set = "Seurat::cc.genes.updated.2019", seurat_version = as.character(packageVersion("Seurat")),
    scoring_normalization = "LogNormalize", scale_factor = 10000, seed = 1, nbin = 24, control_genes_per_marker = min(length(s_genes), length(g2m_genes)), regression_performed = FALSE,
    s_markers_total = length(s_markers), s_markers_matched = length(s_genes), s_markers_detected = s_detected, s_markers_missing = paste(s_missing, collapse = ";"),
    g2m_markers_total = length(g2m_markers), g2m_markers_matched = length(g2m_genes), g2m_markers_detected = g2m_detected, g2m_markers_missing = paste(g2m_missing, collapse = ";"),
    nuclei_G1 = unname(phase_counts["G1"]), nuclei_S = unname(phase_counts["S"]), nuclei_G2M = unname(phase_counts["G2M"]), nuclei_Undecided = unname(phase_counts["Undecided"]),
    percent_G1 = 100 * unname(phase_counts["G1"]) / ncol(object), percent_S = 100 * unname(phase_counts["S"]) / ncol(object), percent_G2M = 100 * unname(phase_counts["G2M"]) / ncol(object), percent_Undecided = 100 * unname(phase_counts["Undecided"]) / ncol(object),
    median_S_score = median(scores$S.Score), median_G2M_score = median(scores$G2M.Score))

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_tsv), recursive = TRUE, showWarnings = FALSE)

saveRDS(object, file = output_rds, compress = FALSE)
write.table(summary, file = summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

message("Scored ", sample_id, ": ", ncol(object), " nuclei; S markers = ", length(s_genes), "/", length(s_markers), "; G2M markers = ", length(g2m_genes), "/", length(g2m_markers), ". Added S.Score, G2M.Score, and Phase only. No regression or filtering.")
