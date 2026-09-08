args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L){stop("Usage: Rscript scripts/05_doublets.R <input.rds> <output.rds> <summary.tsv> <original_combined.h5>", call. = FALSE)}

input_rds <- args[[1]]
output_rds <- args[[2]]
summary_tsv <- args[[3]]
original_h5 <- args[[4]] # use these counts for 'recovered' nuclei (an estimate actually 🤓)

if (!file.exists(input_rds)){stop("Input file doesn't exist: ", input_rds, call. = FALSE)}
if (!file.exists(original_h5)){stop("Original H5 file doesn't exist: ", original_h5, call. = FALSE)}
if (anyDuplicated(normalizePath(c(input_rds, original_h5, output_rds, summary_tsv), mustWork = FALSE))){stop("Input and output paths must be different.", call. = FALSE)}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DoubletFinder))

################ get the count data directly from original h5 

read_original_barcodes <- function(h5_file, sample_id){h5 <- hdf5r::H5File$new(h5_file, mode = "r")
    on.exit(h5$close_all()) 

    if (!h5$exists("matrix/barcodes") || !"library_ids" %in% hdf5r::h5attr_names(h5)){stop("Original H5 requires matrix/barcodes and the library_ids attribute.", call. = FALSE)}
    library_ids <- as.character(hdf5r::h5attr(h5, "library_ids"))
    library_index <- which(library_ids == sample_id)

    if (length(library_index) != 1L){stop("Sample must match 1 original H5 library: ", sample_id, call. = FALSE)}
    original_barcodes <- h5[["matrix/barcodes"]][]
    barcode_groups <- sub("^.*-", "", original_barcodes)

    if (anyNA(original_barcodes) || any(!barcode_groups %in% as.character(seq_along(library_ids)))){stop("Original H5 barcode suffixes do not match its library order.", call. = FALSE)}
    sample_barcodes <- original_barcodes[barcode_groups == as.character(library_index)]
    if (length(sample_barcodes) == 0L || anyDuplicated(sample_barcodes)){stop("Original H5 library must contain unique barcodes: ", sample_id, call. = FALSE)}
    return(sample_barcodes)}

###############################################

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

original_sample_barcodes <- read_original_barcodes(original_h5, sample_id)

if (!all(colnames(object) %in% original_sample_barcodes)){stop("Post-QC barcodes dont all belong to selected original H5 library.", call. = FALSE)}

#H5 barcode count as recovery proxy. I don't have any of the experimental data. how were they loaded, how many recovered, etc (includes all removed by downstream QC)
n_recovered <- length(original_sample_barcodes)
nuclei_post_qc <- ncol(object)
rate_per_1000 <- 0.008 # 0.8%/1,000 recovered nuclei using 10x v2 user guide 
expected_rate <- rate_per_1000 * (n_recovered / 1000)
expected_doublets <- round(expected_rate * nuclei_post_qc)
if (expected_rate >= 1 || expected_doublets < 1L || expected_doublets >= nuclei_post_qc){stop("Expected doublet count must be at least 1 and less than the post-QC nucleus count. Check estimate.", call. = FALSE)}

work_object <- object #copy for preprocessing steps to feed to DoubletFInder

DefaultAssay(work_object) <- "RNA"

set.seed(1234)
work_object <- NormalizeData(work_object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
work_object <- FindVariableFeatures(work_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
work_object <- ScaleData(work_object, features = VariableFeatures(work_object), verbose = FALSE)
work_object <- RunPCA(work_object, features = VariableFeatures(work_object), npcs = 50, seed.use=1234, verbose = FALSE)

### PRINT LOADINGS (8 to 15)
#print(work_object[["pca"]], dims = 8:15, nfeatures = 5) # TO inspect positive and negative loadings

############# START PCA ELBOW and HEATMAPS ####################
###############################################################

#pca_plot <- ElbowPlot(work_object, ndims = 50) + ggplot2::labs(title = paste("PCA elbow:", sample_id))
#ggplot2::ggsave(filename = file.path(dirname(summary_tsv), paste0(sample_id, ".pca_elbow.png")), plot = pca_plot, width = 8, height = 5, dpi = 180)

#pc_heatmaps <- DimHeatmap(work_object, dims = 12:15, nfeatures = 30, cells = 500, reduction = "pca", balanced = TRUE, ncol = 2, fast = FALSE, combine = TRUE)
#ggplot2::ggsave(filename = file.path(dirname(summary_tsv), paste0(sample_id, ".pca_heatmaps_12_15.png")), plot = pc_heatmaps, width = 14, height = 12, units = "in", dpi = 180, bg = "white")

############### END PCA ELBOW and HEATMAPS ####################
###############################################################

pcs <- 1:15

#####################################################################
######################### START Find pK #############################
#####################################################################

dir.create(dirname(summary_tsv), recursive = TRUE, showWarnings = FALSE)

set.seed(1234)
sweep_results <- paramSweep(work_object, PCs = pcs, sct = FALSE, num.cores = 3)


################ START - Comparison of doublet score distribution ##############

# saveRDS(sweep_results, file = file.path(dirname(summary_tsv), paste0(sample_id, ".pK_sweep_results.rds")))

# scores_001 <- sweep_results[["pN_0.25_pK_0.01"]]
# scores_030 <- sweep_results[["pN_0.25_pK_0.3"]]

# if (is.null(scores_001) || is.null(scores_030)){stop("Required pN-pK combinations are missing from sweep.", call. = FALSE)}

# score_data <- rbind(data.frame(Barcode = rownames(scores_001), pK = "0.01", pANN = scores_001$pANN), data.frame(Barcode = rownames(scores_030), pK = "0.30", pANN = scores_030$pANN))

# if (any(!is.finite(score_data$pANN)) || any(score_data$pANN < 0 | score_data$pANN > 1)){stop("pANN scores must be finite values between 0 and 1.", call. = FALSE)}

# write.table(score_data, file = file.path(dirname(summary_tsv), paste0(sample_id, ".pANN_comparison.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

# score_plot <- ggplot2::ggplot(score_data, ggplot2::aes(x = pANN)) +
#     ggplot2::geom_histogram(breaks = seq(0, 1, by = 0.01), fill = "#3C6E91", colour = "white") +
#     ggplot2::facet_wrap(~pK, ncol = 1, labeller = "label_both") +
#     ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
#     ggplot2::theme_bw() +
#     ggplot2::labs(title = paste("pANN distributions:", sample_id), subtitle = "pN = 0.25; PCs = 1:15", x = "pANN score", y = "Number of nuclei")
# ggplot2::ggsave(filename = file.path(dirname(summary_tsv), paste0(sample_id, ".pANN_comparison.png")), plot = score_plot, width = 9, height = 7, dpi = 180, bg = "white")

################ END ##############

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


valid_rows <- is.finite(pk_table$pK) & pk_table$pK > 0 & pk_table$pK < 1 & is.finite(pk_table$BCmetric)

if (!any(valid_rows)){stop("No pK candidates.", call. = FALSE)}
if (any(!valid_rows)){warning("Excluded invalid pK candidates.", call. = FALSE)}

pk_candidates <- pk_table[valid_rows, , drop = FALSE]
peak_rows <- which(pk_candidates$BCmetric == max(pk_candidates$BCmetric))

if (length(peak_rows) != 1L){stop("Multiple pK values share maximum BCmetric. Check sweep!", call. = FALSE)}

best_pK <- pk_candidates$pK[peak_rows]
best_BCmetric <- pk_candidates$BCmetric[peak_rows]
pk_boundary <- best_pK %in% range(pk_candidates$pK)

if (pk_boundary){warning("Selected pK is at boundary of valid sweep range. Check sweep.", call. = FALSE)}

pk_selection <- data.frame(sample_id = sample_id, condition = condition, pcs = paste(pcs, collapse = ","), pK = best_pK, BCmetric = best_BCmetric, boundary = pk_boundary, excluded_candidates = sum(!valid_rows))

write.table(pk_selection, file = file.path(dirname(summary_tsv), paste0(sample_id, ".pK_selected.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Sample:", sample_id, "\nSelected pK:", best_pK, "\nBCmetric:", best_BCmetric, "\n")

###################################################################
######################### END Find pK #############################
###################################################################

########## get counts from original h5s with pK selected from above ############################

expected_summary <- data.frame(sample_id = sample_id, condition = condition, original_h5 = normalizePath(original_h5), count_source = "original_h5_barcodes", n_recovered = n_recovered, nuclei_post_qc = nuclei_post_qc, rate_per_1000 = rate_per_1000, expected_rate = expected_rate, expected_rate_percent = 100 * expected_rate, expected_doublets = expected_doublets, pK = best_pK, BCmetric = best_BCmetric)
write.table(expected_summary, file = file.path(dirname(summary_tsv), paste0(sample_id, ".expected_doublets.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Original H5 barcodes:", n_recovered, "\nPost-QC nuclei:", nuclei_post_qc, "\nExpected rate (%):", 100 * expected_rate, "\nExpected doublets:", expected_doublets, "\nSelected pK:", best_pK, "\n")

#doubletFinder uses pK = best_pK and nExp = expected_doublets ahhhhhhhhhh #######################

set.seed(1234)
work_object <- doubletFinder(work_object, PCs = pcs, pN = 0.25, pK = best_pK, nExp = expected_doublets, reuse.pANN = NULL, sct = FALSE)

classification_column <- paste("DF.classifications", 0.25, best_pK, expected_doublets, sep = "_")

if (!classification_column %in% colnames(work_object[[]])) {stop("DoubletFinder classification column missing.", call. = FALSE)}

singlet_barcodes <- colnames(work_object)[work_object[[]][[classification_column]] == "Singlet"]
work_object <- subset(work_object, cells = singlet_barcodes)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(work_object, file = output_rds)
cat("Doublets removed:", nuclei_post_qc - ncol(work_object), "\n")
cat("Singlets saved:", ncol(work_object), "\n")
cat("Output:", output_rds, "\n")