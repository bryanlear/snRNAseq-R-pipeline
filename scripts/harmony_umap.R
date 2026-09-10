library(Seurat)
library(harmony)
library(ggplot2)

input_path <- "results/GSE174367/umap/brain_pca.rds"
output_dir <- "results/GSE174367/umap"
dims_use <- 1:30

brain_merged <- readRDS(input_path)

DefaultAssay(brain_merged) <- "SCT"
brain_merged$Batch <- factor(brain_merged$Batch)

sample_design <- unique(brain_merged[[]][, c("sample_id", "Batch", "condition"), drop = FALSE])
sample_design <- sample_design[order(sample_design$Batch, sample_design$condition, sample_design$sample_id), , drop = FALSE]

print(sample_design, row.names = FALSE)
print(with(sample_design, table(Batch, condition)))
print(with(brain_merged[[]], table(Batch, condition)))

############## ---- HARMONY COORDINATES ------ ##############

set.seed(12345)
brain_merged <- RunHarmony(object = brain_merged, group.by.vars = "Batch", reduction.use = "pca", dims.use = dims_use, reduction.save = "harmony", project.dim = TRUE, max_iter = 10, early_stop = TRUE, plot_convergence = FALSE, verbose = TRUE)

############## ------ HARMONY UMAP ---------- ##############

brain_merged <- RunUMAP(object = brain_merged, reduction = "harmony", dims = dims_use, n.neighbors = 30L, min.dist = 0.3, reduction.name = "umap.harmony", reduction.key = "UMAPharmony_", seed.use = 12345, verbose = TRUE)
saveRDS(brain_merged, file = file.path(output_dir, "brain_harmony_umap.rds"), compress = FALSE)

#🤬