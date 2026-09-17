args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L){stop("Usage: Rscript scripts/downstream_analysis/histogram_cell_type.R <input.rsd> <output_dir>", call. = FALSE)}

input_rds <- args[[1]]
output_dir <- args[[2]]

library(Seurat)
library(ggplot2)

brain_merged <- readRDS(input_rds)
metadata <- brain_merged[[]]

if (!all(c("sample_id", "Cell.Type") %in% colnames(metadata))){stop("Metadata needs sample_id and Cell.Type.", call. = FALSE)}

# print(table(metadata$sample_id, metadata$Cell.Type, useNA = "ifany")) #sample_id and cell type (Cell.Type in .rds)

cell_type_counts <- table(metadata$sample_id, metadata$Cell.Type)
cell_type_percentages <- 100 * prop.table(cell_type_counts, margin = 1)
print(round(cell_type_percentages, 2))
print(rowSums(cell_type_percentages))

##############----- UMAP COLORS-------################

p_cell_types <- DimPlot(object = brain_merged, reduction = "umap.harmony", group.by = "Cell.Type", label = TRUE, repel = TRUE, pt.size = 1, shuffle = TRUE, seed = 12345, raster = TRUE, raster.dpi = c(1024, 1024), combine = FALSE)[[1]]
umap_color_scale <- ggplot_build(p_cell_types)$plot$scales$get_scales("colour")
cell_type_colors <- setNames(umap_color_scale$map(umap_color_scale$get_limits()), umap_color_scale$get_limits())

stopifnot(all(colnames(cell_type_percentages) %in% names(cell_type_colors)))

############## ----- HISTOGRAM -------################

cell_type_composition <- as.data.frame(as.table(cell_type_percentages), responseName = "percentage")
colnames(cell_type_composition) <- c("sample_id", "cell_type", "percentage")

if (!"Diagnosis" %in% colnames(metadata)){stop("Metadata needs Diagnosis.", call. = FALSE)}

sample_diagnosis <- unique(metadata[, c("sample_id", "Diagnosis"), drop = FALSE])

if (anyNA(sample_diagnosis) || anyDuplicated(sample_diagnosis$sample_id) || !all(sample_diagnosis$Diagnosis %in% c("Control", "AD"))) {stop("Each sample must have one diagnosis: Control or AD.", call. = FALSE)}

cell_type_composition$Diagnosis <- factor(sample_diagnosis$Diagnosis[match(as.character(cell_type_composition$sample_id), sample_diagnosis$sample_id)], levels = c("Control", "AD"))

sample_order <- unique(as.character(cell_type_composition$sample_id))
sample_order <- sample_order[order(as.integer(sub("^Sample-", "", sample_order)))]
cell_type_composition$sample_id <- factor(cell_type_composition$sample_id, levels = sample_order)
cell_type_composition$cell_type <- factor(cell_type_composition$cell_type, levels = names(cell_type_colors))

p_cell_type_composition <- ggplot(cell_type_composition, aes(x = sample_id, y = percentage, fill = cell_type)) +
    geom_col(width = 0.8) +
    facet_grid(cols = vars(Diagnosis), scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cell_type_colors, drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Sample", y = "Cell-type percentage", fill = "Cell type", title = "Cell-type composition per sample") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(filename = file.path(output_dir, "cell_type_composition.png"), plot = p_cell_type_composition, width = max(10, length(sample_order) * 0.3), height = 6, dpi = 300, bg = "white")
