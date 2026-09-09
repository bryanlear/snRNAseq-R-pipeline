library(Seurat)
library(tools)
library(Matrix)
library(ggplot2)

# input_dir <- "results/GSE174367/06_ambient_rna"
# rds_files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))

# if (length(rds_files) == 0L) {stop("No .rds files.")}

# brain_clean <- lapply(rds_files, function(path){object <- readRDS(path)
# stopifnot(inherits(object, "Seurat"), "decontX" %in% Assays(object), all(c("sample_id", "condition", "Batch", "Cell.Type", "cluster") %in% colnames(object[[]])))
#     DefaultAssay(object) <- "decontX"
#     return(object)})

# names(brain_clean) <- file_path_sans_ext(basename(rds_files))
# nuclei_per_sample <- vapply(brain_clean, ncol, numeric(1))
# print(nuclei_per_sample)

# cat("Samples:", length(brain_clean), "\n")
# cat("Total nuclei:", sum(nuclei_per_sample), "\n")

##################--------------------#####################
##################### SCTransform ########################

# brain_sct <- lapply(names(brain_clean), function(sample_id){message("SCTransform: ", sample_id)
# object <- brain_clean[[sample_id]]
#     corrected_totals <- colSums(LayerData(object, assay = "decontX", layer = "counts"))
#     if (any(corrected_totals <= 0)){stop("NO corrected nuclei found in: ", sample_id)}
#     object <- SCTransform(object, assay = "decontX", new.assay.name = "SCT", vst.flavor = "v2", variable.features.n = 3000, vars.to.regress = NULL, conserve.memory = TRUE, seed.use = 12345, verbose = TRUE)
#     return(object)})

# names(brain_sct) <- names(brain_clean)

# sct_summary <- do.call(rbind, lapply(names(brain_sct), function(sample_id){object <- brain_sct[[sample_id]]
# data.frame(sample_id = sample_id, nuclei = ncol(object), variable_features = length(VariableFeatures(object, assay = "SCT")), active_assay = DefaultAssay(object))}))

# print(sct_summary)
# cat("Total nuclei after SCTransform:", sum(sct_summary$nuclei), "\n")

################## ---------- PCA with normalized counts (SCTransform) ----------#####################
################## --------------------------------------------------------------#####################

# pca_features <- SelectIntegrationFeatures(object.list = brain_sct, nfeatures = 3000, verbose = TRUE)
# brain_sct <- PrepSCTIntegration(object.list = brain_sct, assay = "SCT", anchor.features = pca_features, verbose = TRUE)
# brain_merged <- merge(x = brain_sct[[1]], y = brain_sct[-1], add.cell.ids = names(brain_sct), merge.data = TRUE, merge.dr = FALSE, project = "GSE174367")

# DefaultAssay(brain_merged) <- "SCT"
# VariableFeatures(brain_merged) <- pca_features

# brain_merged <- RunPCA(object = brain_merged, assay = "SCT", features = pca_features, npcs = 50, reduction.name = "pca", seed.use = 12345, verbose = TRUE)

# cat("Selected genes:", length(pca_features), "\n")
# cat("Merged nuclei:", ncol(brain_merged), "\n")
# cat("PCA assay:", DefaultAssay(brain_merged[["pca"]]), "\n")
# print(dim(Embeddings(brain_merged, reduction = "pca")))

# output_dir <- "results/GSE174367/umap"
# dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ggsave(filename = file.path(output_dir, "pca_elbow.png"), plot = ElbowPlot(brain_merged, reduction = "pca", ndims = 50), width = 8, height = 5, dpi = 180)
# saveRDS(brain_merged, file = file.path(output_dir, "brain_pca.rds"), compress = FALSE)

################## ---------------------- First UMAP ----------------------------#####################
################## --------------------------------------------------------------#####################

output_dir <- "results/GSE174367/umap"
brain_merged <- readRDS(file.path(output_dir, "brain_pca.rds"))

dims_use <- 1:30

brain_merged <- RunUMAP(object = brain_merged, reduction = "pca", dims = dims_use, n.neighbors = 30L, min.dist = 0.3, reduction.name = "umap.unintegrated", reduction.key = "UMAPraw_", seed.use = 12345, verbose = TRUE)

print(dim(Embeddings(brain_merged, reduction = "umap.unintegrated")))

saveRDS(brain_merged, file = file.path(output_dir, "brain_umap_unintegrated.rds"), compress = FALSE)

p_cell_types <- DimPlot(object = brain_merged, reduction = "umap.unintegrated", group.by = "Cell.Type", label = TRUE, repel = TRUE, pt.size = 1, shuffle = TRUE, seed = 12345, raster = TRUE, raster.dpi = c(1024, 1024)) + ggtitle("Cell types before batch correction")
p_structure <- DimPlot(object = brain_merged, reduction = "umap.unintegrated", group.by = c("Cell.Type", "sample_id", "Batch", "condition"), ncol = 2, pt.size = 1, shuffle = TRUE, seed = 12345, raster = TRUE, raster.dpi = c(1024, 1024))

ggsave(filename = file.path(output_dir, "umap_cell_types.png"), plot = p_cell_types, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(filename = file.path(output_dir, "umap_structure.png"), plot = p_structure, width = 16, height = 12, dpi = 300, bg = "white")
