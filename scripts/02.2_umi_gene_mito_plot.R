args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L) {stop("Usage: Rscript scripts/02.2_umi_gene_mito_plot.R <cells.tsv> <umi_mito.png> <gene_mito.png> <umi_gene_mito.png>", call. = FALSE)}

cells_tsv <- args[[1]]
output_files <- args[2:4]

if (!file.exists(cells_tsv)) {stop("Input file doesn't exist, homie: ", cells_tsv, call. = FALSE)}
if (anyDuplicated(output_files)) {stop("Use different output path for each plot.", call. = FALSE)}

suppressPackageStartupMessages(library(ggplot2))

cells <- read.delim(cells_tsv, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("sample_id", "nCount_RNA", "nFeature_RNA", "percent.mt")

if (nrow(cells) == 0L || !all(required %in% names(cells))) {stop("Input requires nuclei and columns: sample_id, nCount_RNA, nFeature_RNA, percent.mt.", call. = FALSE)}

sample_id <- unique(as.character(cells$sample_id))

if (length(sample_id) != 1L || anyNA(sample_id) || !nzchar(sample_id)) {stop("Input must contain one nonempty sample_id.", call. = FALSE)}
if (!all(vapply(cells[, required[-1], drop = FALSE], function(x) is.numeric(x) && all(is.finite(x)), logical(1)))) {stop("Counts and mitochondrial percentages must be finite numeric values.", call. = FALSE)}
if (any(cells$nCount_RNA <= 0 | cells$nFeature_RNA <= 0)) {stop("Counts must be positive for logarithmic axes.", call. = FALSE)}
if (any(cells$percent.mt < 0 | cells$percent.mt > 100)) {stop("Mitochondrial percentages must be between 0 and 100.", call. = FALSE)}

#Nuclei with higher mitochondrial percentages are last so they're visible
cells <- cells[order(cells$percent.mt), , drop = FALSE]
subtitle <- paste0(sample_id, " | ", format(nrow(cells), big.mark = ","), " nuclei | one point per nucleus")
plot_theme <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0))
umi_mito <- ggplot(cells, aes(x = nCount_RNA, y = percent.mt)) + geom_point(colour = "#245B82", size = 1, alpha = 0.6) + scale_x_log10(labels = scales::label_comma()) + scale_y_continuous(limits = c(0, NA)) + labs(title = "UMI counts vs mitochondrial RNA", subtitle = subtitle, x = "Total UMIs per nucleus (log10 scale)", y = "Mitochondrial RNA (%)", caption = "All input nuclei are shown. No filtering threshold is applied.") + plot_theme
gene_mito <- ggplot(cells, aes(x = nFeature_RNA, y = percent.mt)) + geom_point(colour = "#245B82", size = 1, alpha = 0.6) + scale_x_log10(labels = scales::label_comma()) + scale_y_continuous(limits = c(0, NA)) + labs(title = "Detected genes vs mitochondrial RNA", subtitle = subtitle, x = "Detected genes per nucleus (log10 scale)", y = "Mitochondrial RNA (%)", caption = "All input nuclei are shown. No filtering threshold is applied.") + plot_theme
umi_gene_mito <- ggplot(cells, aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) + geom_point(size = 1.2, alpha = 0.8) + scale_x_log10(labels = scales::label_comma()) + scale_y_log10(labels = scales::label_comma()) + scale_colour_viridis_c(option = "plasma", name = "Mitochondrial\nRNA (%)") + labs(title = "UMIs, detected genes, and mitochondrial RNA", subtitle = subtitle, x = "Total UMIs per nucleus (log10 scale)", y = "Detected genes per nucleus (log10 scale)", caption = "All input nuclei are shown. Color limits are specific to this sample.") + plot_theme

plots <- list(umi_mito, gene_mito, umi_gene_mito)

for (i in seq_along(output_files)) {
  dir.create(dirname(output_files[[i]]), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = output_files[[i]], plot = plots[[i]], device = "png", width = 8, height = 5.5, units = "in", dpi = 180, bg = "white")}

message("Plotted ", sample_id, ": ", nrow(cells), " nuclei --> 3 PNG files. No nuclei removed.")
