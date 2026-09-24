# Rscript scripts/geneformer_downstream_analysis/knn.R input.csv output_dir [batch_size]
args <- commandArgs(trailingOnly = TRUE)

if (!length(args) %in% c(2L, 3L)) stop("Usage: Rscript scripts/geneformer_downstream_analysis/knn.R <input.csv> <output_dir> [batch_size]", call. = FALSE)

input_csv <- args[[1L]]
output_dir <- args[[2L]]
batch_size <- if (length(args) == 3L) suppressWarnings(as.integer(args[[3L]])) else 100L # default batches of 100; batches of 1000 take about 16 minutes for 53,000 nuclei
# The Euclidean pass takes about 17 additional minutes with the same data.

embedding_dimension <- 768L
index_column_number <- 1L
embedding_column_numbers <- seq.int(2L, embedding_dimension + 1L)
first_metadata_column_number <- embedding_dimension + 2L
k_values <- c(15L, 30L, 50L)
distance_metrics <- c("cosine", "euclidean")
required_packages <- c("data.table", "FNN", "igraph")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]

if (length(missing_packages) > 0L) stop(sprintf("Required package(s) needed: %s", paste(missing_packages, collapse = ", ")), call. = FALSE)
if (!file.exists(input_csv)) stop(sprintf("Input file doesn't exist: %s", input_csv), call. = FALSE)
if (is.na(batch_size) || batch_size < 1L) stop("batch_size must be a positive integer.", call. = FALSE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Read input CSV: ", input_csv)
input_table <- data.table::fread(input_csv, check.names = FALSE)

if (ncol(input_table) < first_metadata_column_number) stop(sprintf("Input has %d columns. It must contain an index column, 768 embedding columns, and at least one metadata column.", ncol(input_table)), call. = FALSE)
if (nrow(input_table) <= max(k_values)) stop(sprintf("Input has %d nuclei. It must have more than the largest k value (%d).", nrow(input_table), max(k_values)), call. = FALSE)

metadata_column_numbers <- seq.int(first_metadata_column_number, ncol(input_table))
embedding_table <- input_table[, ..embedding_column_numbers]
metadata <- as.data.frame(input_table[, ..metadata_column_numbers], check.names = FALSE)
input_index <- input_table[[index_column_number]]
required_metadata_columns <- c("sample_id", "condition", "batch", "cell_type")
missing_metadata_columns <- setdiff(required_metadata_columns, names(metadata))
numeric_embedding_columns <- vapply(embedding_table, is.numeric, logical(1L))

if (length(missing_metadata_columns) > 0L) stop(sprintf("Required metadata column(s) missing after embedding column 767: %s", paste(missing_metadata_columns, collapse = ", ")), call. = FALSE)
if (!all(numeric_embedding_columns)) stop(sprintf("Embedding columns 2:769 must be numeric. Non-numeric column(s): %s", paste(names(embedding_table)[!numeric_embedding_columns], collapse = ", ")), call. = FALSE)

embedding_matrix <- as.matrix(embedding_table)
storage.mode(embedding_matrix) <- "double"

if (anyNA(embedding_matrix) || any(!is.finite(embedding_matrix))) stop("Embedding columns contain NA, NaN, or infinite values.", call. = FALSE)

row_id <- seq_len(nrow(input_table))
metadata_output <- data.frame(.knn_row_id = row_id, .input_index = input_index, metadata, check.names = FALSE)
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
    graph <- igraph::set_vertex_attr(graph, ".input_index", value = input_index)
    graph <- igraph::set_graph_attr(graph, "distance_metric", metric)
    graph <- igraph::set_graph_attr(graph, "k", as.integer(k))
    graph <- igraph::set_graph_attr(graph, "directed_knn", TRUE)
    graph <- igraph::set_graph_attr(graph, "embedding_columns", names(embedding_table))
    graph <- igraph::set_graph_attr(graph, "vertex_metadata", metadata_output)
    graph}
run_batched_knn <- function(search_matrix, metric, maximum_k, batch_size) {
    nucleus_count <- nrow(search_matrix)
    neighbor_index <- matrix(NA_integer_, nrow = nucleus_count, ncol = maximum_k)
    neighbor_distance <- matrix(NA_real_, nrow = nucleus_count, ncol = maximum_k)
    batch_starts <- seq.int(1L, nucleus_count, by = batch_size)
    search_started_at <- Sys.time()
    message(sprintf("Run %s KNN search: %d nuclei, maximum k = %d, batch size = %d", metric, nucleus_count, maximum_k, batch_size))

    for (batch_start in batch_starts) {
        batch_end <- min(batch_start + batch_size - 1L, nucleus_count)
        batch_rows <- seq.int(batch_start, batch_end)
        batch_result <- FNN::get.knnx(data = search_matrix, query = search_matrix[batch_rows, , drop = FALSE], k = maximum_k + 1L, algorithm = "brute")

        for (batch_row in seq_along(batch_rows)) {
            nucleus_row <- batch_rows[[batch_row]]
            nonself_positions <- which(batch_result$nn.index[batch_row, ] != nucleus_row)
            if (length(nonself_positions) < maximum_k) stop(sprintf("Could not find %d non-self neighbors for nucleus row %d.", maximum_k, nucleus_row), call. = FALSE)
            selected_positions <- nonself_positions[seq_len(maximum_k)]
            neighbor_index[nucleus_row, ] <- batch_result$nn.index[batch_row, selected_positions]
            neighbor_distance[nucleus_row, ] <- batch_result$nn.dist[batch_row, selected_positions]}

        completed <- batch_end
        elapsed_seconds <- as.numeric(difftime(Sys.time(), search_started_at, units = "secs"))
        remaining_seconds <- if (completed < nucleus_count) (elapsed_seconds / completed) * (nucleus_count - completed) else 0
        message(sprintf("[%s] %d/%d nuclei (%.1f%%) | elapsed %.1f min | estimated remaining %.1f min", metric, completed, nucleus_count, 100 * completed / nucleus_count, elapsed_seconds / 60, remaining_seconds / 60))}

    list(nn.index = neighbor_index, nn.dist = neighbor_distance)}

normalized_embedding_matrix <- normalize_rows(embedding_matrix)
manifest_rows <- vector("list", length(distance_metrics) * length(k_values))
manifest_position <- 1L
input_file_information <- file.info(input_csv)
input_signature <- list(path = normalizePath(input_csv), size = unname(input_file_information$size), modified = as.character(input_file_information$mtime))

for (metric in distance_metrics) {
    search_matrix <- if (metric == "cosine") normalized_embedding_matrix else embedding_matrix
    search_cache_file <- file.path(output_dir, sprintf("knn_%s_k%d_search_cache.rds", metric, max(k_values)))
    cached_search <- if (file.exists(search_cache_file)) readRDS(search_cache_file) else NULL
    cache_is_valid <- !is.null(cached_search) && identical(cached_search$input_signature, input_signature) && identical(cached_search$metric, metric) && identical(cached_search$maximum_k, max(k_values)) && identical(dim(cached_search$result$nn.index), c(nrow(input_table), max(k_values)))
    if (cache_is_valid) message("Resume from completed search cache: ", search_cache_file)
    if (!cache_is_valid) {
        knn_result <- run_batched_knn(search_matrix = search_matrix, metric = metric, maximum_k = max(k_values), batch_size = batch_size)
        saveRDS(list(input_signature = input_signature, metric = metric, maximum_k = max(k_values), result = knn_result), search_cache_file, compress = FALSE)
        message("Write completed search cache: ", search_cache_file)} else knn_result <- cached_search$result
    all_neighbor_index <- knn_result$nn.index
    all_neighbor_distance <- knn_result$nn.dist

    if (metric == "cosine") {
        all_neighbor_distance <- (all_neighbor_distance * all_neighbor_distance) / 2
        all_neighbor_distance[all_neighbor_distance < 0] <- 0
        all_neighbor_distance[all_neighbor_distance > 2] <- 2}

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
