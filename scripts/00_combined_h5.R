args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L){stop("Usage: Rscript scripts/00_combined_h5.R <input.h5> <metadata.csv.gz> <output_dir> <error|exclude>", call. = FALSE)}

input_h5 <- args[[1]]
metadata_file <- args[[2]]
output_dir <- args[[3]]
unmatched_policy <- args[[4]] #error to stop if any [barcode] lacks metadata or exlude to get rid of those empty barcodes (I use this one)

if (!unmatched_policy %in% c("error", "exclude")){stop("Unmatched policy must be error or exclude.", call. = FALSE)}

suppressPackageStartupMessages(library(Seurat))

metadata <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("Barcode", "SampleID", "Diagnosis")

if (!all(required %in% names(metadata))){stop("Metadata requires Barcode, SampleID, and Diagnosis.", call. = FALSE)}
if (nrow(metadata) == 0L || anyNA(metadata[required]) || any(trimws(as.matrix(metadata[required])) == "")){stop("Required metadata values cannot be empty or missing.", call. = FALSE)}
if (anyDuplicated(metadata$Barcode)){stop("Metadata barcodes must be unique.", call. = FALSE)}

sample_ids <- unique(metadata$SampleID)

if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", sample_ids))){stop("Sample IDs contain characters unsuitable for filenames.", call. = FALSE)}

conditions_per_sample <- tapply(metadata$Diagnosis, metadata$SampleID, function(x) length(unique(x)))

if (any(conditions_per_sample != 1L)){stop("Each sample must have 1 diagnosis.", call. = FALSE)}

counts <- Read10X_h5(input_h5, use.names = TRUE)

if (is.list(counts)){stop("Expected one gene-expression matrix.", call. = FALSE)}

barcodes <- colnames(counts)

if (anyDuplicated(barcodes)) {stop("H5 barcodes must be unique.", call. = FALSE)}
if (any(!metadata$Barcode %in% barcodes)){stop("Some metadata barcodes are missing from the H5.", call. = FALSE)}

unmatched <- setdiff(barcodes, metadata$Barcode)

message(length(barcodes), " H5 barcodes; ", nrow(metadata), " metadata rows; ", length(unmatched), " barcodes without metadata.")

if (length(unmatched) > 0L && unmatched_policy == "error"){stop("H5 barcodes lack metadata. Use exclude only if exclusion is intended.", call. = FALSE)}

#metadata must match H5 barcode order.
matched_barcodes <- barcodes[barcodes %in% metadata$Barcode]
metadata <- metadata[match(matched_barcodes, metadata$Barcode), , drop = FALSE]
rownames(metadata) <- metadata$Barcode
#so they column names match downstream steps
metadata$sample_id <- metadata$SampleID
metadata$condition <- metadata$Diagnosis

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.table(data.frame(Barcode = unmatched, reason = rep("missing_metadata", length(unmatched))), file = file.path(output_dir, "excluded_barcodes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sample_ids <- unique(metadata$sample_id)
summaries <- vector("list", length(sample_ids))

for (i in seq_along(sample_ids)){sample_id <- sample_ids[[i]]
    cells <- rownames(metadata)[metadata$sample_id == sample_id]

    object <- CreateSeuratObject(counts = counts[, cells, drop = FALSE], project = sample_id, min.cells = 0, min.features = 0, meta.data = metadata[cells, , drop = FALSE])

    saveRDS(object, file = file.path(output_dir, paste0(sample_id, ".rds")), compress = FALSE)

    summaries[[i]] <- data.frame(sample_id = sample_id, condition = unique(object$condition), nuclei = ncol(object))

    message("Saved ", sample_id, ": ", ncol(object), " nuclei.")
    rm(object)}

write.table(do.call(rbind, summaries), file = file.path(output_dir, "import_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)