################################################################################
# Illustrative smoothclust and boundary-density analyses for de Oliveira P2CRC
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(irlba)
library(Matrix)
library(SpatialExperiment)
library(smoothclust)

dataset_id <- "deOliveira"
processed_dir <- file.path("data", dataset_id, "processed")
results_dir <- file.path("results", dataset_id)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 123
n_threads <- 8
smoothing_k <- 18
n_pcs <- 15
sketch_n <- 30000

settings <- list(
  VisiumHD = list(
    file = "spe_deOliveira_VisiumHD_P2CRC_8um.rds",
    platform = "Visium HD",
    spatial_unit = "8-micron square bin",
    n_clusters = 15,
    boundary_k = 8
  ),
  Xenium = list(
    file = "spe_deOliveira_Xenium_P2CRC.rds",
    platform = "Xenium",
    spatial_unit = "segmented cell",
    n_clusters = 15,
    boundary_k = 6
  )
)

elapsed <- function(start) {
  as.numeric(difftime(Sys.time(), start, units = "secs"))
}

for (platform_id in names(settings)) {
  cfg <- settings[[platform_id]]
  total_start <- Sys.time()

  stage_start <- Sys.time()
  spe <- readRDS(file.path(processed_dir, cfg$file))
  read_seconds <- elapsed(stage_start)
  coords <- spatialCoords(spe)
  input_nnz <- length(assay(spe, "counts")@x)

  stage_start <- Sys.time()
  # Normalize with full-library totals before smoothing and log transformation.
  spe <- smoothclust(
    spe,
    method = "knn",
    k = smoothing_k,
    n_threads = n_threads,
    size_factors = colData(spe)$library_size
  )
  smoothing_seconds <- elapsed(stage_start)
  smoothed <- assay(spe, "logcounts_smooth")
  smoothed_nnz <- length(smoothed@x)

  set.seed(seed)
  sketch <- sample(seq_len(ncol(spe)), min(sketch_n, ncol(spe)))
  stage_start <- Sys.time()
  # Fit PCA on the sketch, then project all locations with the same centering.
  pca_fit <- prcomp_irlba(
    t(smoothed[, sketch, drop = FALSE]),
    n = n_pcs,
    center = TRUE,
    scale. = FALSE
  )
  pca <- as.matrix(t(smoothed) %*% pca_fit$rotation)
  pca_center <- as.numeric(pca_fit$center %*% pca_fit$rotation)
  pca <- sweep(pca, 2, pca_center, "-")
  pca_seconds <- elapsed(stage_start)

  stage_start <- Sys.time()
  set.seed(seed)
  km <- kmeans(
    pca[sketch, , drop = FALSE],
    centers = cfg$n_clusters,
    nstart = 3,
    iter.max = 500,
    algorithm = "Lloyd"
  )
  # Maximizing these scores minimizes squared distance to the fitted centroids.
  scores <- 2 * tcrossprod(pca, km$centers)
  scores <- sweep(scores, 2, rowSums(km$centers ^ 2), "-")
  smooth_labels <- factor(
    max.col(scores, ties.method = "first"),
    levels = seq_len(cfg$n_clusters)
  )
  clustering_seconds <- elapsed(stage_start)
  rm(scores, pca, pca_fit, km)

  stage_start <- Sys.time()
  smooth_bd <- boundary_density(
    coords,
    smooth_labels,
    k = cfg$boundary_k,
    n_threads = n_threads,
    adjust = "analytic"
  )
  smooth_bd_seconds <- elapsed(stage_start)

  author_bd <- NULL
  author_bd_seconds <- NA_real_
  if (platform_id == "VisiumHD") {
    stage_start <- Sys.time()
    author_bd <- boundary_density(
      coords,
      colData(spe)$author_cluster,
      k = cfg$boundary_k,
      n_threads = n_threads,
      adjust = "analytic"
    )
    author_bd_seconds <- elapsed(stage_start)
  }

  spot_results <- data.frame(
    location_id = colnames(spe),
    x = coords[, 1],
    y = coords[, 2],
    smoothclust_label = smooth_labels,
    local_bd_smoothclust = smooth_bd$local_boundary_density
  )
  if (platform_id == "VisiumHD") {
    spot_results$author_cluster <- colData(spe)$author_cluster
    spot_results$author_deconvolution <- colData(spe)$decon_label_1
    spot_results$local_bd_author_cluster <-
      author_bd$local_boundary_density
  }
  saveRDS(
    spot_results,
    file.path(
      results_dir,
      paste0("deOliveira_", platform_id, "_spot_results.rds")
    )
  )

  bd_summary <- data.frame(
    dataset_id = dataset_id,
    platform = cfg$platform,
    sample_id = "P2CRC",
    spatial_unit = cfg$spatial_unit,
    label_source = c(
      "smoothclust + sketched PCA/k-means",
      if (platform_id == "VisiumHD") "OSTA distributed clusters"
    ),
    n_locations = ncol(spe),
    n_genes = nrow(spe),
    n_labels = c(
      length(unique(smooth_labels)),
      if (platform_id == "VisiumHD") {
        length(unique(colData(spe)$author_cluster))
      }
    ),
    smoothing_method = c(
      "knn",
      if (platform_id == "VisiumHD") NA_character_
    ),
    smoothing_k = c(
      smoothing_k,
      if (platform_id == "VisiumHD") NA_integer_
    ),
    boundary_k = cfg$boundary_k,
    boundary_density = c(
      smooth_bd$boundary_density,
      if (platform_id == "VisiumHD") author_bd$boundary_density
    ),
    expected_boundary_density = c(
      smooth_bd$expected_boundary_density,
      if (platform_id == "VisiumHD") {
        author_bd$expected_boundary_density
      }
    ),
    relative_boundary_density = c(
      smooth_bd$relative_boundary_density,
      if (platform_id == "VisiumHD") {
        author_bd$relative_boundary_density
      }
    )
  )
  write.csv(
    bd_summary,
    file.path(
      results_dir,
      paste0("deOliveira_", platform_id, "_BD_summary.csv")
    ),
    row.names = FALSE
  )

  runtime <- data.frame(
    dataset_id = dataset_id,
    platform = cfg$platform,
    sample_id = "P2CRC",
    n_locations = ncol(spe),
    n_genes = nrow(spe),
    input_nnz = input_nnz,
    smoothed_nnz = smoothed_nnz,
    n_threads = n_threads,
    stage = c(
      "read_processed_object",
      "smoothclust",
      "sketched_pca_and_full_projection",
      "sketched_kmeans_and_full_assignment",
      "boundary_density_smoothclust",
      if (platform_id == "VisiumHD") "boundary_density_author_clusters"
    ),
    seconds = c(
      read_seconds,
      smoothing_seconds,
      pca_seconds,
      clustering_seconds,
      smooth_bd_seconds,
      if (platform_id == "VisiumHD") author_bd_seconds
    )
  )

  nn_distances <- BiocNeighbors::findKNN(
    coords,
    k = max(smoothing_k, cfg$boundary_k),
    get.index = FALSE,
    get.distance = "transposed",
    num.threads = n_threads
  )$distance
  # These thresholds flag long edges; they do not remove graph neighbors.
  distance_threshold <- if (platform_id == "VisiumHD") 5 else 50
  coordinate_unit <- if (platform_id == "VisiumHD") {
    "8-micron bin-coordinate units"
  } else {
    "micrometers"
  }
  graph_k <- c(boundary_density = cfg$boundary_k, smoothing = smoothing_k)
  distance_summary <- do.call(
    rbind,
    lapply(
      names(graph_k),
      function(graph) {
        k <- unname(graph_k[[graph]])
        kth_distance <- nn_distances[k, ]
        data.frame(
          graph = graph,
          k = k,
          coordinate_unit = coordinate_unit,
          median_kth_neighbor_distance = median(kth_distance),
          q99_kth_neighbor_distance = unname(quantile(kth_distance, 0.99)),
          max_kth_neighbor_distance = max(kth_distance),
          long_edge_threshold = distance_threshold,
          n_locations_above_threshold = sum(kth_distance > distance_threshold)
        )
      }
    )
  )
  write.csv(
    distance_summary,
    file.path(
      results_dir,
      paste0("deOliveira_", platform_id, "_neighbor_distance_summary.csv")
    ),
    row.names = FALSE
  )

  total_row <- runtime[1, , drop = FALSE]
  total_row$stage <- "total_including_neighbor_diagnostics"
  total_row$seconds <- elapsed(total_start)
  runtime <- rbind(runtime, total_row)
  write.csv(
    runtime,
    file.path(
      results_dir,
      paste0("deOliveira_", platform_id, "_runtime.csv")
    ),
    row.names = FALSE
  )

  rm(spe, smoothed, spot_results, smooth_bd, author_bd)
}
