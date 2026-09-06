args <- commandArgs(trailingOnly = TRUE)

if (!length(args) %in% c(6L, 7L)){stop("Usage: Rscript scripts/03_mit_contamination.R <input.rds> <output.rds> <summary.tsv> <fixed|mad> <gene_pattern> <cutoff_percent|nmads> [min_diff]", call. = FALSE)}

#fixed '^MT-' 5 percent
#mad '^MT-' 3 0.5 (3 * mad_mt, 0.5) 

input_rds <- args[[1]]
output_rds <- args[[2]]
summary_tsv <- args[[3]]

method <- args[[4]]
mito_pattern <- args[[5]] #-^MT
value <- suppressWarnings(as.numeric(args[[6]]))
min_diff <- if (length(args) == 7L) suppressWarnings(as.numeric(args[[7]])) else NA_real_

if (!method %in% c("fixed", "mad")){stop("Method must be fixed or mad.", call. = FALSE)}
if (method == "fixed" && length(args) != 6L){stop("Fixed mode needs six arguments.", call. = FALSE)}
if (method == "mad" && length(args) != 7L){stop("MAD mode needs seven arguments including min_diff.", call. = FALSE)}
if (!is.finite(value)){stop("The cutoff or nmads must be a finite number.", call. = FALSE)}
if (method == "fixed" && (value < 0 || value > 100)) {stop("Fixed cutoff must be between 0 and 100 percent.", call. = FALSE)}
if (method == "mad" && value <= 0) {stop("nmads must be positive.", call. = FALSE)}
if (method == "mad" && (!is.finite(min_diff) || min_diff < 0 || min_diff > 100)){stop("min_diff must be between 0 and 100 percentage points.", call. = FALSE)}
if (!nzchar(mito_pattern)){stop("Gene pattern can't be empty.", call. = FALSE)}
if (!file.exists(input_rds)){stop("Input file doesnt exist: ", input_rds, call. = FALSE)}

suppressPackageStartupMessages(library(Seurat))

object <- readRDS(input_rds)

if (!inherits(object, "Seurat") || ncol(object) == 0L){stop("Input must be a Seurat object containing nuclei.", call. = FALSE)}
if (!"RNA" %in% Assays(object)) {stop("The RNA assay is required.", call. = FALSE)}
if (!all(c("sample_id", "condition") %in% colnames(object[[]]))){stop("Metadata requires sample_id and condition.", call. = FALSE)}

sample_id <- unique(as.character(object$sample_id))
condition <- unique(as.character(object$condition))

if (length(sample_id) != 1L || anyNA(sample_id) || !nzchar(sample_id)){stop("Input must contain one nonempty sample_id.", call. = FALSE)}
if (length(condition) != 1L || anyNA(condition) || !nzchar(condition)){stop("Input must contain one nonempty condition.", call. = FALSE)}

mito_genes <- grep(mito_pattern, rownames(object[["RNA"]]), value = TRUE)

if (length(mito_genes) == 0L){stop("No mitochondrial genes match the supplied pattern.", call. = FALSE)}

object$percent.mt <- PercentageFeatureSet(object, features = mito_genes, assay = "RNA")
percent_mt <- object$percent.mt

if (anyNA(percent_mt) || any(!is.finite(percent_mt)) || any(percent_mt < 0 | percent_mt > 100)) {stop("Mitochondrial percentages must be finite values between 0 and 100.", call. = FALSE)}

median_mt <- median(percent_mt)
mad_mt <- mad(percent_mt)

if (method == "fixed") {threshold <- value} else {threshold <- median_mt + max(value * mad_mt, min_diff)}

fail_mt <- percent_mt > threshold
keep <- !fail_mt

message("Sample: ", sample_id, "; method: ", method, "; threshold: ", signif(threshold, 6), "%; flagged: ", sum(fail_mt), " of ", length(keep), ".")

if (!any(keep)){stop("The selected threshold removes every nucleus. Review the setting before proceeding.", call. = FALSE)}

filtered_object <- subset(object, cells = colnames(object)[keep])

summary <- data.frame(
    sample_id = sample_id, condition = condition, method = method, mito_pattern = mito_pattern,
    mitochondrial_features = length(mito_genes), fixed_cutoff_percent = if (method == "fixed") value else NA_real_,
    nmads = if (method == "mad") value else NA_real_, min_diff_percentage_points = min_diff,
    median_percent_mt = median_mt, mad_percent_mt = mad_mt, threshold_percent_mt = threshold,
    nuclei_before = length(keep), nuclei_removed = sum(fail_mt), nuclei_after = sum(keep),
    percent_removed = round(100 * mean(fail_mt), 2))

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_tsv), recursive = TRUE, showWarnings = FALSE)

saveRDS(filtered_object, file = output_rds, compress = FALSE)
write.table(summary, file = summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

message("Saved ", sample_id, ": ", sum(keep), " nuclei retained.")