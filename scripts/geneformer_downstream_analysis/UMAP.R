# Rscript scripts/geneformer_downstream_analysis/UMAP.R results/geneformer_downstream_analysis/knn_output results/geneformer_downstream_analysis/umap_output
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) stop("Usage: Rscript scripts/geneformer_downstream_analysis/UMAP.R <knn_input_dir> <umap_output_dir>", call. = FALSE)

input_dir <- args[[1L]]
output_dir <- args[[2L]]

distance_metrics <- c("cosine", "euclidean")
k_values <- c(15L, 30L, 50L)
random_seed <- 42L
minimum_distance <- 0.1
required_packages <- c("data.table", "ggplot2", "ggrepel", "igraph", "uwot")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]

if (length(missing_packages) > 0L) stop(sprintf("Required package(s) needed: %s", paste(missing_packages, collapse = ", ")), call. = FALSE)
if (!dir.exists(input_dir)) stop(sprintf("KNN input directory does not exist: %s", input_dir), call. = FALSE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expected_graph_files <- unlist(lapply(distance_metrics, function(metric) file.path(input_dir, sprintf("knn_%s_k%d.rds", metric, k_values))), use.names = FALSE)
missing_graph_files <- expected_graph_files[!file.exists(expected_graph_files)]

if (length(missing_graph_files) > 0L) stop(sprintf("Required KNN graph file(s) missing: %s", paste(basename(missing_graph_files), collapse = ", ")), call. = FALSE)

extract_precomputed_neighbors <- function(graph, expected_k) {
    nucleus_count <- igraph::vcount(graph)
    edge_table <- data.table::as.data.table(igraph::as_edgelist(graph, names = FALSE))
    data.table::setnames(edge_table, c("from_row_id", "to_row_id"))
    edge_table[, neighbor_rank := as.integer(igraph::edge_attr(graph, "neighbor_rank"))]
    edge_table[, distance := as.numeric(igraph::edge_attr(graph, "distance"))]
    if (igraph::ecount(graph) != nucleus_count * expected_k) stop(sprintf("Graph has %d edges; expected %d nuclei x %d neighbors = %d edges.", igraph::ecount(graph), nucleus_count, expected_k, nucleus_count * expected_k), call. = FALSE)
    if (anyNA(edge_table$distance) || any(!is.finite(edge_table$distance)) || any(edge_table$distance < 0)) stop("Graph contains invalid neighbor distances.", call. = FALSE)
    if (anyDuplicated(edge_table[, .(from_row_id, neighbor_rank)])) stop("Graph contains duplicate neighbor ranks for a nucleus.", call. = FALSE)
    if (!setequal(unique(edge_table$neighbor_rank), seq_len(expected_k))) stop(sprintf("Graph neighbor ranks must contain 1:%d.", expected_k), call. = FALSE)
    data.table::setorder(edge_table, from_row_id, neighbor_rank)
    neighbor_index <- matrix(edge_table$to_row_id, nrow = nucleus_count, ncol = expected_k, byrow = TRUE)
    neighbor_distance <- matrix(edge_table$distance, nrow = nucleus_count, ncol = expected_k, byrow = TRUE)
    umap_index <- cbind(seq_len(nucleus_count), neighbor_index)
    umap_distance <- cbind(rep(0, nucleus_count), neighbor_distance)
    list(idx = umap_index, dist = umap_distance)}

make_umap_plot <- function(plot_table, color_column, title, facet_by_condition = FALSE) {
    plot_data <- data.table::copy(plot_table)
    plot_data[, (color_column) := factor(as.character(get(color_column)))]
    group_levels <- levels(droplevels(plot_data[[color_column]]))
    group_colors <- stats::setNames(grDevices::hcl.colors(length(group_levels), palette = "Dark 3"), group_levels)
    label_group_columns <- if (facet_by_condition) unique(c("condition", color_column)) else color_column
    label_table <- plot_data[!is.na(get(color_column)), .(UMAP_1 = stats::median(UMAP_1), UMAP_2 = stats::median(UMAP_2)), by = label_group_columns]
    legend_title <- gsub("_", " ", color_column, fixed = TRUE)
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = .data[[color_column]])) + ggplot2::geom_point(size = 0.35, alpha = 0.75, stroke = 0) + ggrepel::geom_label_repel(data = label_table, ggplot2::aes(x = UMAP_1, y = UMAP_2, label = .data[[color_column]], color = .data[[color_column]]), inherit.aes = FALSE, fill = "white", fontface = "bold", size = 3.5, box.padding = 0.4, point.padding = 0.2, min.segment.length = 0, max.overlaps = Inf, seed = random_seed, show.legend = FALSE) + ggplot2::scale_color_manual(values = group_colors, na.value = "grey70", drop = FALSE) + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4, alpha = 1))) + ggplot2::labs(title = title, color = legend_title, x = "UMAP 1", y = "UMAP 2") + ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "right", legend.key.width = grid::unit(0.8, "cm"), legend.key.height = grid::unit(0.6, "cm"), plot.title = ggplot2::element_text(face = "bold"))
    if (facet_by_condition) plot <- plot + ggplot2::facet_wrap(ggplot2::vars(condition))
    plot}

manifest_rows <- vector("list", length(expected_graph_files))
manifest_position <- 1L

for (metric in distance_metrics) {
    for (k in k_values) {
        graph_file <- file.path(input_dir, sprintf("knn_%s_k%d.rds", metric, k))
        message(sprintf("Read KNN graph: %s", graph_file))
        graph <- readRDS(graph_file)

        if (!inherits(graph, "igraph")) stop(sprintf("File is not an igraph object: %s", graph_file), call. = FALSE)
        if (!identical(as.character(igraph::graph_attr(graph, "distance_metric")), metric)) stop(sprintf("Graph metric does not match filename: %s", graph_file), call. = FALSE)
        if (!identical(as.integer(igraph::graph_attr(graph, "k")), k)) stop(sprintf("Graph k does not match filename: %s", graph_file), call. = FALSE)
        metadata <- igraph::graph_attr(graph, "vertex_metadata")
        if (!is.data.frame(metadata) || nrow(metadata) != igraph::vcount(graph)) stop(sprintf("Graph metadata is missing or has the wrong row count: %s", graph_file), call. = FALSE)
        missing_metadata_columns <- setdiff(c("condition", "cell_type"), names(metadata))
        if (length(missing_metadata_columns) > 0L) stop(sprintf("Graph metadata column(s) missing: %s", paste(missing_metadata_columns, collapse = ", ")), call. = FALSE)

        file_prefix <- sprintf("umap_%s_k%d", metric, k)
        umap_file <- file.path(output_dir, paste0(file_prefix, ".rds"))
        cached_umap <- if (file.exists(umap_file)) readRDS(umap_file) else NULL
        cache_is_valid <- !is.null(cached_umap) && identical(cached_umap$metric, metric) && identical(as.integer(cached_umap$k), k) && identical(as.integer(cached_umap$seed), random_seed) && identical(as.numeric(cached_umap$min_dist), minimum_distance) && identical(dim(cached_umap$coordinates), c(igraph::vcount(graph), 2L))
        precomputed_neighbors <- NULL

        if (cache_is_valid) {message(sprintf("Reuse existing UMAP coordinates for %s k=%d", metric, k))
            umap_coordinates <- cached_umap$coordinates} else {
            precomputed_neighbors <- extract_precomputed_neighbors(graph, k)
            message(sprintf("Run UMAP for %s k=%d with %d nuclei and %d stored neighbors", metric, k, igraph::vcount(graph), k))
            set.seed(random_seed)
            umap_coordinates <- uwot::umap(X = NULL, n_components = 2L, min_dist = minimum_distance, init = "spectral", nn_method = precomputed_neighbors, n_sgd_threads = 1L, verbose = TRUE, seed = random_seed)}
        plot_table <- data.table::as.data.table(data.table::copy(metadata))
        plot_table[, UMAP_1 := umap_coordinates[, 1L]]
        plot_table[, UMAP_2 := umap_coordinates[, 2L]]
        plot_table[, condition_cell_type := interaction(condition, cell_type, sep = " | ", drop = TRUE)]
        coordinate_file <- file.path(output_dir, paste0(file_prefix, "_coordinates.csv"))
        condition_plot_file <- file.path(output_dir, paste0(file_prefix, "_condition.png"))
        cell_type_plot_file <- file.path(output_dir, paste0(file_prefix, "_cell_type.png"))
        condition_cell_type_plot_file <- file.path(output_dir, paste0(file_prefix, "_condition_cell_type.png"))
        data.table::fwrite(plot_table, coordinate_file)
        saveRDS(list(coordinates = umap_coordinates, metadata = metadata, metric = metric, k = k, seed = random_seed, min_dist = minimum_distance), umap_file, compress = FALSE)
        condition_plot <- make_umap_plot(plot_table, "condition", sprintf("UMAP: %s distance, k = %d, colored by condition", metric, k))
        cell_type_plot <- make_umap_plot(plot_table, "cell_type", sprintf("UMAP: %s distance, k = %d, colored by cell type", metric, k))
        condition_cell_type_plot <- make_umap_plot(plot_table, "cell_type", sprintf("UMAP: %s distance, k = %d, cell type by condition", metric, k), facet_by_condition = TRUE)
        ggplot2::ggsave(condition_plot_file, condition_plot, width = 10, height = 8, units = "in", dpi = 300, bg = "white")
        ggplot2::ggsave(cell_type_plot_file, cell_type_plot, width = 10, height = 8, units = "in", dpi = 300, bg = "white")
        ggplot2::ggsave(condition_cell_type_plot_file, condition_cell_type_plot, width = 14, height = 9, units = "in", dpi = 300, bg = "white")
        manifest_rows[[manifest_position]] <- data.frame(distance_metric = metric, k = k, nuclei = nrow(plot_table), coordinate_file = basename(coordinate_file), umap_file = basename(umap_file), condition_plot = basename(condition_plot_file), cell_type_plot = basename(cell_type_plot_file), condition_cell_type_plot = basename(condition_cell_type_plot_file), seed = random_seed, min_dist = minimum_distance, stringsAsFactors = FALSE)
        manifest_position <- manifest_position + 1L
        message(sprintf("Completed UMAP for %s k=%d", metric, k))
        rm(graph, metadata, precomputed_neighbors, umap_coordinates, plot_table, condition_plot, cell_type_plot, condition_cell_type_plot)
        invisible(gc())}}

manifest <- data.table::rbindlist(manifest_rows)
data.table::fwrite(manifest, file.path(output_dir, "manifest.csv"))
message("All UMAP outputs are complete: ", normalizePath(output_dir))
