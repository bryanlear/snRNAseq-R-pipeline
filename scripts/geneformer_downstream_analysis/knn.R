# Rscript scripts/geneformer_downstream_analysis/knn.R input.csv output_dir
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) stop("Usage: Rscript scripts/geneformer_downstream_analysis/knn.R <input.csv> <output_dir>", call. = FALSE)

input_csv <- args[[1L]]
output_dir <- args[[2L]]


embedding_dimension <- 768L
k_values <- c(15L, 30L, 50L)
distance_metrics <- c("cosine", "euclidean")
required_packages <- c("data.table", "FNN", "igraph")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]

if (length(missing_packages) > 0L) stop(sprintf("Required package(s) needed: %s", paste(missing_packages, collapse = ", ")), call. = FALSE)
if (!file.exists(input_csv)) stop(sprintf("Input file doesn't exist: %s", input_csv), call. = FALSE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Read input CSV: ", input_csv)
input_table <- data.table::fread(input_csv, check.names = FALSE)

if (ncol(input_table) < embedding_dimension) stop(sprintf("Input has %d columns. It must have at least %d embedding columns.", ncol(input_table), embedding_dimension), call. = FALSE)
if (nrow(input_table) <= max(k_values)) stop(sprintf("Input has %d nuclei. It must have more than the largest k value (%d).", nrow(input_table), max(k_values)), call. = FALSE)

embedding_column_numbers <- seq_len(embedding_dimension)
metadata_column_numbers <- if (ncol(input_table) > embedding_dimension) seq.int(embedding_dimension + 1L, ncol(input_table)) else integer(0L)
embedding_table <- input_table[, ..embedding_column_numbers]
numeric_embedding_columns <- vapply(embedding_table, is.numeric, logical(1L))

if (!all(numeric_embedding_columns)) stop(sprintf("Embedding columns 2:768 must be numeric. Non-numeric column(s): %s", paste(names(embedding_table)[!numeric_embedding_columns], collapse = ", ")), call. = FALSE)

embedding_matrix <- as.matrix(embedding_table)
storage.mode(embedding_matrix) <- "double"

if (anyNA(embedding_matrix) || any(!is.finite(embedding_matrix))) stop("Embedding columns contain NA, NaN, or infinite values.", call. = FALSE)

metadata <- if (length(metadata_column_numbers) > 0L) as.data.frame(input_table[, ..metadata_column_numbers], check.names = FALSE) else data.frame(row.names = seq_len(nrow(input_table)))
row_id <- seq_len(nrow(input_table))
metadata_output <- data.frame(.knn_row_id = row_id, metadata, check.names = FALSE)
data.table::fwrite(metadata_output, file.path(output_dir, "metadata.csv"))
normalize_rows <- function(matrix) {row_norms <- sqrt(rowSums(matrix * matrix))
    if (any(row_norms == 0)) stop("Cosine distance is undefined for an embedding with zero L2 norm.", call. = FALSE)
    matrix / row_norms}
make_graph <- function(neighbor_index, neighbor_distance, metric, k) {
    edge_table <- data.table::data.table(from_row_id = rep(row_id, each = k), to_row_id = as.vector(t(neighbor_index)), neighbor_rank = rep(seq_len(k), times = length(row_id)), distance = as.vector(t(neighbor_distance)))
    edge_vector <- as.vector(t(as.matrix(edge_table[, .(from_row_id, to_row_id)])))
    graph <- igraph::make_empty_graph(n = nrow(input_table), directed = TRUE)
    graph <- igraph::add_edges(graph, edge_vector)
    graph <- igraph::set_edge_attr(graph, "distance", value = edge_table$distance)
    graph <- igraph::set_edge_attr(graph, "neighbor_rank", value = edge_table$neighbor_rank)
    graph <- igraph::set_vertex_attr(graph, ".knn_row_id", value = row_id)
    graph <- igraph::set_graph_attr(graph, "distance_metric", metric)
    graph <- igraph::set_graph_attr(graph, "k", as.integer(k))
    graph <- igraph::set_graph_attr(graph, "directed_knn", TRUE)
    graph <- igraph::set_graph_attr(graph, "embedding_columns", names(embedding_table))
    graph <- igraph::set_graph_attr(graph, "vertex_metadata", metadata_output)
    graph}

normalized_embedding_matrix <- normalize_rows(embedding_matrix)
manifest_rows <- vector("list", length(distance_metrics) * length(k_values))
manifest_position <- 1L

for (metric in distance_metrics) {
    search_matrix <- if (metric == "cosine") normalized_embedding_matrix else embedding_matrix
    message("Run ", metric, " KNN search with maximum k = ", max(k_values))
    knn_result <- FNN::get.knn(data = search_matrix, k = max(k_values), algorithm = "brute") #then i may try k-d tree, NN descent, balltree
    all_neighbor_index <- knn_result$nn.index
    all_neighbor_distance <- knn_result$nn.dist

    if (metric == "cosine") all_neighbor_distance <- pmax(0, pmin(2, (all_neighbor_distance * all_neighbor_distance) / 2))

    for (k in k_values) {
        neighbor_index <- all_neighbor_index[, seq_len(k), drop = FALSE]
        neighbor_distance <- all_neighbor_distance[, seq_len(k), drop = FALSE]
        graph <- make_graph(neighbor_index, neighbor_distance, metric, k)
        graph_file <- file.path(output_dir, sprintf("knn_%s_k%d.rds", metric, k))
        saveRDS(graph, graph_file, compress = FALSE)
        manifest_rows[[manifest_position]] <- data.frame(distance_metric = metric, k = k, graph_file = basename(graph_file), nuclei = nrow(input_table), dimensions = embedding_dimension, directed_edges = nrow(input_table) * k, stringsAsFactors = FALSE)
        manifest_position <- manifest_position + 1L
        message("Write graph: ", graph_file)}}

manifest <- data.table::rbindlist(manifest_rows)
data.table::fwrite(manifest, file.path(output_dir, "manifest.csv"))
message("KNN graph done: ", normalizePath(output_dir))

#😩