args <- commandArgs(trailingOnly = TRUE)

if (!length(args) %in% 5:7) {stop("Usage: Rscript scripts/02.2_umi_gene_mito_plot.R <cells.tsv> <umi_mito.png> <gene_mito.png> <umi_gene_mito.png> <flagged.tsv> [nmads=3] [min_diff=0.5]", call. = FALSE)}

cells_tsv <- args[[1]]
output_files <- args[2:4]
flagged_tsv <- args[[5]]
nmads <- if (length(args) >= 6L) suppressWarnings(as.numeric(args[[6]])) else 3
min_diff <- if (length(args) == 7L) suppressWarnings(as.numeric(args[[7]])) else 0.5

if (!file.exists(cells_tsv)) {stop("Input file doesn't exist, homie: ", cells_tsv, call. = FALSE)}
if (anyDuplicated(normalizePath(c(cells_tsv, output_files, flagged_tsv), mustWork = FALSE))) {stop("Use different paths for the input and every output.", call. = FALSE)}
if (!is.finite(nmads) || nmads <= 0) {stop("nmads must be a positive, finite number.", call. = FALSE)}
if (!is.finite(min_diff) || min_diff < 0 || min_diff > 100) {stop("min_diff must be between 0 and 100 percentage points.", call. = FALSE)}

suppressPackageStartupMessages(library(ggplot2))

cells <- read.delim(cells_tsv, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("sample_id", "nCount_RNA", "nFeature_RNA", "percent.mt", "Barcode")

if (nrow(cells) == 0L || !all(required %in% names(cells))) {stop("Input requires nuclei and columns: Barcode, sample_id, nCount_RNA, nFeature_RNA, percent.mt.", call. = FALSE)}
if (anyNA(cells$Barcode) || any(!nzchar(trimws(cells$Barcode))) || anyDuplicated(cells$Barcode)) {stop("Barcodes must be nonempty and unique.", call. = FALSE)}

sample_id <- unique(as.character(cells$sample_id))

if (length(sample_id) != 1L || anyNA(sample_id) || !nzchar(sample_id)) {stop("Input must contain one nonempty sample_id.", call. = FALSE)}
if (!all(vapply(cells[, c("nCount_RNA", "nFeature_RNA", "percent.mt"), drop = FALSE], function(x) is.numeric(x) && all(is.finite(x)), logical(1)))) {stop("Counts and mitochondrial percentages must be finite numeric values.", call. = FALSE)}
if (any(cells$nCount_RNA <= 0 | cells$nFeature_RNA <= 0)) {stop("Counts must be positive for logarithmic axes.", call. = FALSE)}
if (any(cells$percent.mt < 0 | cells$percent.mt > 100)) {stop("Mitochondrial percentages must be between 0 and 100.", call. = FALSE)}

median_mt <- median(cells$percent.mt)
mad_mt <- mad(cells$percent.mt)
raw_mad_threshold <- median_mt + nmads * mad_mt
threshold <- median_mt + max(nmads * mad_mt, min_diff)

if (!is.finite(threshold)) {stop("The calculated threshold is not finite. Check nmads.", call. = FALSE)}

cells$high_mt <- cells$percent.mt > threshold
flagged <- cells[cells$high_mt, , drop = FALSE]
flagged$nmads <- rep(nmads, nrow(flagged))
flagged$min_diff_percentage_points <- rep(min_diff, nrow(flagged))
flagged$median_percent_mt <- rep(median_mt, nrow(flagged))
flagged$mad_percent_mt <- rep(mad_mt, nrow(flagged))
flagged$raw_mad_threshold_percent_mt <- rep(raw_mad_threshold, nrow(flagged))
flagged$threshold_percent_mt <- rep(threshold, nrow(flagged))
flagged <- flagged[order(flagged$percent.mt, decreasing = TRUE), , drop = FALSE]

#Nuclei with higher mitochondrial percentages are last so they're visible
cells <- cells[order(cells$percent.mt), , drop = FALSE]
subtitle <- paste0(sample_id, " | ", format(nrow(cells), big.mark = ","), " nuclei | one point per nucleus\nUpper cutoff = ", signif(threshold, 5), "% | flagged: ", nrow(flagged), " (", round(100 * nrow(flagged) / nrow(cells), 2), "%)")
threshold_caption <- paste0("Cutoff = median + max(", nmads, " x MAD, ", min_diff, " percentage points). R MAD scaling: 1.4826.")
plot_theme <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0))
umi_mito <- ggplot(cells, aes(x = nCount_RNA, y = percent.mt)) + geom_point(colour = "#245B82", size = 1, alpha = 0.6) + geom_hline(yintercept = threshold, colour = "#B32424", linetype = "dashed", linewidth = 0.6) + scale_x_log10(labels = scales::label_comma()) + scale_y_continuous(limits = c(0, NA)) + labs(title = "UMI counts vs mitochondrial RNA", subtitle = subtitle, x = "Total UMIs per nucleus (log10 scale)", y = "Mitochondrial RNA (%)", caption = paste0(threshold_caption, "\nDashed line: upper cutoff. Values above it are flagged. No nuclei removed.")) + plot_theme
gene_mito <- ggplot(cells, aes(x = nFeature_RNA, y = percent.mt)) + geom_point(colour = "#245B82", size = 1, alpha = 0.6) + geom_hline(yintercept = threshold, colour = "#B32424", linetype = "dashed", linewidth = 0.6) + scale_x_log10(labels = scales::label_comma()) + scale_y_continuous(limits = c(0, NA)) + labs(title = "Detected genes vs mitochondrial RNA", subtitle = subtitle, x = "Detected genes per nucleus (log10 scale)", y = "Mitochondrial RNA (%)", caption = paste0(threshold_caption, "\nDashed line: upper cutoff. Values above it are flagged. No nuclei removed.")) + plot_theme
umi_gene_mito <- ggplot(cells, aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) + geom_point(size = 1.2, alpha = 0.8) + geom_point(data = flagged, shape = 21, fill = NA, colour = "black", size = 2.2, stroke = 0.5) + scale_x_log10(labels = scales::label_comma()) + scale_y_log10(labels = scales::label_comma()) + scale_colour_viridis_c(option = "plasma", name = "Mitochondrial\nRNA (%)") + labs(title = "UMIs, detected genes, and mitochondrial RNA", subtitle = subtitle, x = "Total UMIs per nucleus (log10 scale)", y = "Detected genes per nucleus (log10 scale)", caption = paste0(threshold_caption, "\nBlack outlines: flagged nuclei. Color limits are sample-specific. No nuclei removed.")) + plot_theme

plots <- list(umi_mito, gene_mito, umi_gene_mito)

for (i in seq_along(output_files)) {
  dir.create(dirname(output_files[[i]]), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = output_files[[i]], plot = plots[[i]], device = "png", width = 8, height = 5.5, units = "in", dpi = 180, bg = "white")}

dir.create(dirname(flagged_tsv), recursive = TRUE, showWarnings = FALSE)
write.table(flagged, file = flagged_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

message("Plotted ", sample_id, ": ", nrow(cells), " nuclei; threshold = ", signif(threshold, 6), "%; flagged = ", nrow(flagged), "; three PNG files and one flagged TSV. No nuclei removed.")
