args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L){stop(paste(
      "Usage: Rscript scripts/01_create_seurat.R",
      "<sample_id> <condition> <input.h5> <output.rds>"),
    call. = FALSE)}

sample_id <- args[[1]]
condition <- args[[2]]
input_h5 <- args[[3]]
output_rds <- args[[4]]

if (!file.exists(input_h5)){stop("Input file doesn't exist: ", input_h5, call. = FALSE)}
if (!condition %in% c("control", "treated")){stop("Unexpected condition: ", condition, call. = FALSE)}

suppressPackageStartupMessages(library(Seurat))
message("Current sample ", sample_id, " from ", input_h5)

counts <- Read10X_h5(filename = input_h5, use.names = TRUE)

if (is.list(counts)){stop("Read10X_h5 returned multiple matrices - expected 1 gene expression matrix.",call. = FALSE)}

object <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 0, min.features = 0)

object$sample_id <- sample_id
object$condition <- condition

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(object,file = output_rds, compress = FALSE)

message("Saved ", sample_id, ": ", nrow(object), " features × ", ncol(object), " nuclei to ",output_rds)
