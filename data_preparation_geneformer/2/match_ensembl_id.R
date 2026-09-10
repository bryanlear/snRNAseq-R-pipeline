library(Seurat)
library(hdf5r)

input_dir <- "results/GSE174367/06_ambient_rna"
output_dir <- "data_preparation_geneformer/geneformer_rds"
h5_file <- "brain_data_human/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"

rds_files <- sort(list.files(input_dir, pattern = "^Sample-[0-9]+\\.rds$", full.names = TRUE))
stopifnot(length(rds_files) == 18L)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

########### ----- FEATURES DATA FROM ORIGINAL ------ ##############
###################################################################

#h5dump -n 1 .h5
#h5dump -d "matrix/features/x" .h5 | head -n 30

h5 <- H5File$new(h5_file, mode = "r")

ensembl_versioned <- h5[["matrix/features/id"]][]
gene_symbol <- h5[["matrix/features/name"]][]

library_ids <- as.character(h5attr(h5, "library_ids"))

h5$close_all()

expected_seurat_names <- make.unique(gene_symbol)
ensembl_id <- sub("_PAR_Y$", "", sub("\\.[0-9]+", "", ensembl_versioned))
seurat_ensembl_names <- gsub("_", "-", ensembl_versioned, fixed = TRUE)

stopifnot(length(ensembl_versioned) == 58721L, !anyDuplicated(ensembl_versioned), !anyDuplicated(seurat_ensembl_names), all(grepl("^ENSG[0-9]+$", ensembl_id)))

feature_metadata <- data.frame(ensembl_id = ensembl_id, ensembl_id_versioned = ensembl_versioned, gene_symbol = gene_symbol, row.names = seurat_ensembl_names, check.names = FALSE)


##################################################################################
########### ----- APPEND ORIGINAL FEATURES TO PROCESSED DATA ------ ##############

for (input_file in rds_files){sample_name <- sub("\\.rds$", "", basename(input_file))

    output_file <- file.path(output_dir, basename(input_file))
    
    if (file.exists(output_file)){stop("Output already exists: ", output_file)} 
    message("Processing ", sample_name)

    object <- readRDS(input_file)

    corrected_counts <- LayerData(object, assay = "decontX", layer = "counts")

    sample_id <- unique(as.character(object$sample_id))
    library_index <- match(sample_id, library_ids)
    barcode_suffix <- sub("^.*-", "", colnames(corrected_counts))

    stopifnot(length(sample_id) == 1L, !is.na(library_index), identical(rownames(corrected_counts), expected_seurat_names), identical(colnames(corrected_counts), rownames(object[[]])), all(barcode_suffix == as.character(library_index)), all(corrected_counts@x >= 0), all(corrected_counts@x == round(corrected_counts@x)))

    restored_counts <- corrected_counts
    rownames(restored_counts) <- seurat_ensembl_names

    stopifnot(identical(restored_counts@x, corrected_counts@x), identical(colnames(restored_counts), colnames(corrected_counts)))

    restored_assay <- CreateAssay5Object(counts = restored_counts)
    restored_assay <- AddMetaData(restored_assay, metadata = feature_metadata)

    object[["decontX"]] <- restored_assay
    DefaultAssay(object) <- "decontX"

    final_counts <- LayerData(object, assay = "decontX", layer = "counts")

    stopifnot(identical(dim(final_counts), dim(corrected_counts)), identical(rownames(final_counts), seurat_ensembl_names), identical(colnames(final_counts), colnames(corrected_counts)), identical(final_counts@x, corrected_counts@x))

    saveRDS(object, output_file, compress = FALSE)
    message("Saved ", sample_name, ": ", nrow(final_counts), " features x ", ncol(final_counts), " nuclei")
    rm(object, corrected_counts, restored_counts, restored_assay, final_counts)
    gc()}

#fuck this shit
