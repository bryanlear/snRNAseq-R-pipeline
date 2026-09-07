args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L){stop("Usage: Rscript scripts/05_doublets.R <input.rds> <output.rds> <summary.tsv>", call. = FALSE)}

input_rds <- args[[1]]
output_rds <- args[[2]]
summary_tsv <- args[[3]]

if (!file.exists(input_rds)){stop("Input file doesn't exist: ", input_rds, call. = FALSE)}
if (anyDuplicated(normalizePath(c(input_rds, output_rds, summary_tsv), mustWork = FALSE))){stop("Input and output paths must be different.", call. = FALSE)}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DoubletFinder))

object <- readRDS(input_rds)

if (!inherits(object, "Seurat")){stop("Input must be Seurat object.", call. = FALSE)}
if (ncol(object) == 0L){stop("Input contains no nuclei.", call. = FALSE)}
if (!"RNA" %in% Assays(object)){stop("Input needs an RNA assay.", call. = FALSE)}

metadata <- object[[]]

if (!all(c("sample_id", "condition") %in% names(metadata))){stop("Metadata requires sample_id and condition.", call. = FALSE)}

sample_id <- unique(as.character(metadata$sample_id))
condition <- unique(as.character(metadata$condition))

if (length(sample_id) != 1L || anyNA(sample_id) || !nzchar(trimws(sample_id))){stop("Input must contain 1 nonempty sample_id.", call. = FALSE)}
if (length(condition) != 1L || anyNA(condition) || !nzchar(trimws(condition))){stop("Input must contain 1 nonempty condition.", call. = FALSE)}

work_object <- object #copy for preprocessing steps to feed to DoubletFInder

DefaultAssay(work_object) <- "RNA"

set.seed(1234)
work_object <- NormalizeData(work_object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
work_object <- FindVariableFeatures(work_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
work_object <- ScaleData(work_object, features = VariableFeatures(work_object), verbose = FALSE)
work_object <- RunPCA(work_object, features = VariableFeatures(work_object), npcs = 50, seed.use=1234, verbose = FALSE)

### PRINT LOADINGS (8 to 15)
print(work_object[["pca"]], dims = 8:15, nfeatures = 5) # TO inspect positive and negative loadings

############# START PCA ELBOW and HEATMAPS ####################
###############################################################

#pca_plot <- ElbowPlot(work_object, ndims = 50) + ggplot2::labs(title = paste("PCA elbow:", sample_id))
#dir.create(dirname(summary_tsv), recursive = TRUE, showWarnings = FALSE)
#ggplot2::ggsave(filename = file.path(dirname(summary_tsv), paste0(sample_id, ".pca_elbow.png")), plot = pca_plot, width = 8, height = 5, dpi = 180)

#pc_heatmaps <- DimHeatmap(work_object, dims = 12:15, nfeatures = 30, cells = 500, reduction = "pca", balanced = TRUE, ncol = 2, fast = FALSE, combine = TRUE)
#ggplot2::ggsave(filename = file.path(dirname(summary_tsv), paste0(sample_id, ".pca_heatmaps_12_15.png")), plot = pc_heatmaps, width = 14, height = 12, units = "in", dpi = 180, bg = "white")

############### END PCA ELBOW and HEATMAPS ####################
###############################################################

pcs <- 1:15

######################### START Find pK #############################

set.seed(1234)
sweep_results <- paramSweep(work_object, PCs = pcs, sct = FALSE, num.cores = 1)
sweep_stats <- summarizeSweep(sweep_results, GT = FALSE)
grDevices::pdf(file = NULL)

pk_table <- find.pK(sweep_stats)

grDevices::dev.off()

pk_table$pK <- as.numeric(as.character(pk_table$pK))

write.table(pk_table, file = file.path(dirname(summary_tsv), paste0(sample_id, ".pK_sweep.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

pk_plot <- ggplot2::ggplot(pk_table, ggplot2::aes(x = pK, y = BCmetric)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = paste("DoubletFinder pK sweep:", sample_id), x = "pK", y = "BCmetric")
ggplot2::ggsave(filename = file.path(dirname(summary_tsv), paste0(sample_id, ".pK_sweep.png")), plot = pk_plot, width = 8, height = 5, dpi = 180, bg = "white")

######################### END Find pK #############################