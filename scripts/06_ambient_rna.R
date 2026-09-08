args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) {stop("Usage: Rscript scripts/06_ambient_rna.R <input.rds> <output.rds>", call. = FALSE)}

input_rds <- args[[1]]
output_rds <- args[[2]]

seed <- 12345L
if (!file.exists(input_rds)){stop("Input file doesnt exist: ", input_rds, call. = FALSE)}
if (anyDuplicated(normalizePath(c(input_rds, output_rds), mustWork = FALSE))) {stop("Input and output paths must be different.", call. = FALSE)}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(decontX))

object <- readRDS(input_rds)

if (!inherits(object, "Seurat") || ncol(object) == 0L){stop("Input must be Seurat object with nuclei.", call. = FALSE)}
if (!"RNA" %in% Assays(object)) {stop("Input requires an RNA assay.", call. = FALSE)}
if ("decontX" %in% Assays(object)) {stop("A decontX assay already exists. Use the output from step 05.", call. = FALSE)}
if (!"sample_id" %in% names(object[[]])) {stop("Metadata requires sample_id.", call. = FALSE)}

sample_id <- unique(as.character(object$sample_id))

if (length(sample_id) != 1L || anyNA(sample_id) || !nzchar(trimws(sample_id))) {stop("Input must contain 1 nonempty sample_id.", call. = FALSE)}

#############################################################
######################## RNA counts ########################
############################################################

counts_layers <- Layers(object[["RNA"]], search = "^counts($|\\.)")

if (!identical(counts_layers, "counts")){stop("RNA requires one counts layer named counts.", call. = FALSE)}

original_counts <- LayerData(object, assay = "RNA", layer = "counts")

if (!inherits(original_counts, "dgCMatrix")){stop("RNA counts must be a dgCMatrix.", call. = FALSE)}
if (!identical(colnames(original_counts), colnames(object))) {stop("RNA counts and object barcodes must have the same order.", call. = FALSE)}
if (any(!is.finite(original_counts@x)) || any(original_counts@x < 0) || any(original_counts@x != round(original_counts@x))) {stop("RNA counts must contain finite, nonnegative integers.", call. = FALSE)}
if (any(Matrix::colSums(original_counts) <= 0)){stop("Each input nucleus must have a positive total count.", call. = FALSE)}

############################################################
######################## DecontX ########################
############################################################

decontx_result <- decontX::decontX(original_counts, z = NULL, background = NULL, seed = seed, verbose = TRUE)

if (!identical(dimnames(decontx_result$decontXcounts), dimnames(original_counts))) {stop("Corrected counts must preserve gene and barcode order.", call. = FALSE)}

contamination <- decontx_result$contamination

if (length(contamination) != ncol(object) || any(!is.finite(contamination)) || any(contamination < 0 | contamination > 1)) {stop("Contamination estimates must contain 1 finite fraction between 0 and 1 per nucleus.", call. = FALSE)}

corrected_counts <- Matrix::drop0(round(decontx_result$decontXcounts))

if (any(!is.finite(corrected_counts@x)) || any(corrected_counts@x < 0)){stop("Corrected counts must be finite and nonnegative.", call. = FALSE)}

######################################################################
######################## Add corrected assay ########################
######################################################################

object[["decontX"]] <- CreateAssay5Object(counts = corrected_counts)

object$decontX_contamination <- setNames(100 * contamination, colnames(original_counts))
object$nCount_decontX <- Matrix::colSums(corrected_counts)
object$nFeature_decontX <- Matrix::colSums(corrected_counts > 0)
DefaultAssay(object) <- "decontX"

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(object, file = output_rds, compress = FALSE)

message("Ambient RNA correction complete: ", sample_id, "; ", ncol(object), " nuclei saved.")