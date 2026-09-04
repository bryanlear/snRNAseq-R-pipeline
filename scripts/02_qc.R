args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L){stop("Usage: Rscript scripts/02_qc.R", "<input.rds> <output.rds> <summary.tsv> <nmads>", call. = FALSE)}

input_rds <- args[[1]]
output_rds <- args[[2]]
summary_tsv <- args[[3]]
nmads <- suppressWarnings(as.numeric(args[[4]])) #number median absolute deviations (MADs) as fourth argument

if (!file.exists(input_rds)){stop("File doesn't exist: ", input_rds, call. = FALSE)}
if (length(nmads) != 1L || !is.finite(nmads) || nmads <= 0){stop("nmads must be positive", call. = FALSE)}

suppressPackageStartupMessages(library(Seurat))

object <- readRDS(input_rds)
metadata <- object[[]]

required_columns <- c("sample_id", "condition", "nCount_RNA", "nFeature_RNA")
missing_columns <- setdiff(required_columns, colnames(metadata))

if (length(missing_columns) > 0L){stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)}

sample_ids <- unique(as.character(object$sample_id))
conditions <- unique(as.character(object$condition))

if (length(sample_ids) != 1L || is.na(sample_ids)){stop("Object must contain one sample_id.", call. = FALSE)}
if (length(conditions) != 1L || is.na(conditions)){stop("Object must contain one condition.", call. = FALSE)}

calculate_limits <- function(values, nmads){if (anyNA(values) || any(!is.finite(values)) || any(values <= 0)){stop("QC values must be positive finite numbers.", call. = FALSE)}

  log_values <- log10(values)
  center <- median(log_values)
  spread <- mad(log_values)

  if (!is.finite(spread) || spread == 0){stop("CANNOT calculate MAD limits because MAD = zero.", call. = FALSE)}

  c(lower = ceiling(10^(center - nmads * spread)), upper = floor(10^(center + nmads * spread)))}

counts <- object$nCount_RNA
features <- object$nFeature_RNA

count_limits <- calculate_limits(counts, nmads)
feature_limits <- calculate_limits(features, nmads)

fail_low_count <- counts < count_limits[["lower"]]
fail_high_count <- counts > count_limits[["upper"]]
fail_low_feature <- features < feature_limits[["lower"]]
fail_high_feature <- features > feature_limits[["upper"]]

keep <- !(
  fail_low_count |
  fail_high_count |
  fail_low_feature |
  fail_high_feature)

filtered_object <- subset(object,cells = colnames(object)[keep])

summary <- data.frame(
  sample_id = sample_ids,
  condition = conditions,
  nmads = nmads,
  count_lower = count_limits[["lower"]],
  count_upper = count_limits[["upper"]],
  feature_lower = feature_limits[["lower"]],
  feature_upper = feature_limits[["upper"]],
  fail_low_count = sum(fail_low_count),
  fail_high_count = sum(fail_high_count),
  fail_low_feature = sum(fail_low_feature),
  fail_high_feature = sum(fail_high_feature),
  nuclei_before = length(keep),
  nuclei_removed = sum(!keep),
  nuclei_after = sum(keep),
  percent_removed = round(100 * sum(!keep) / length(keep), 2))

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_tsv), recursive = TRUE, showWarnings = FALSE)
saveRDS(filtered_object, file = output_rds, compress = FALSE)
write.table(summary, file = summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("QC complete", sample_ids, ": ", length(keep), " nuclei before; ", sum(keep), " after; ", sum(!keep), " removed.")
