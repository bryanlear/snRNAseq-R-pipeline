library(Seurat)
library(ggplot2)

output_dir <- "results/GSE174367/umap"
brain_merged <- readRDS(file.path(output_dir, "brain_harmony_umap.rds"))

p_harmony_cell_types <- DimPlot(object = brain_merged, reduction = "umap.harmony", group.by = "Cell.Type", label = TRUE, repel = TRUE, pt.size = 1, shuffle = TRUE, seed = 12345, raster = TRUE, raster.dpi = c(1024, 1024)) + ggtitle("Cell types after Harmony batch correction")
p_harmony_structure <- DimPlot(object = brain_merged, reduction = "umap.harmony", group.by = c("Cell.Type", "sample_id", "Batch", "condition"), ncol = 2, pt.size = 1, shuffle = TRUE, seed = 12345, raster = TRUE, raster.dpi = c(1024, 1024))

ggsave(filename = file.path(output_dir, "umap_harmony_cell_types.png"), plot = p_harmony_cell_types, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(filename = file.path(output_dir, "umap_harmony_structure.png"), plot = p_harmony_structure, width = 16, height = 12, dpi = 300, bg = "white")