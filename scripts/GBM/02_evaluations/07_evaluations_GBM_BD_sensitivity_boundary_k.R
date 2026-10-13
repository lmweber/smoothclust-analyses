################################################################################
# Script to evaluate BD neighborhood size sensitivity for boundary density
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(smoothclust)
library(SpatialExperiment)

dataset_id <- "GBM"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)

bd_input_file <- file.path(results_dir, "GBM_BD_summary.csv")
spe_input_file <- file.path(results_dir, "spe_GBM_BD_plot.rds")
spe_plot_file <- file.path(
  results_dir,
  "spe_GBM_BD_sensitivity_boundary_k_plot.rds"
)

smoothing_method <- "uniform"
smoothing_bandwidth <- 0.05
boundary_k_values <- c(2L, 4L, 6L, 8L, 12L, 18L, 24L, 36L)
default_boundary_k <- 6L
cluster_numbers <- 3:8
parallel_threads <- 14

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# helper functions

make_label_col <- function(n_clusters) {
  paste0("smoothclust_k", n_clusters)
}

make_local_col <- function(label_col, boundary_k) {
  paste0("local_bd_", label_col, "_bdk", sprintf("%02d", boundary_k))
}


# load data

bd_input <- read.csv(bd_input_file)
spe <- readRDS(spe_input_file)

sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))


# recompute BD summaries while holding smoothclust labels fixed

spe_list <- list()
bd_summary <- data.frame()

for (sample_id in sample_ids) {
  message("Processing sample: ", sample_id)
  
  spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
  spe_out <- spe_sub
  geo_sample_id <- unique(as.character(colData(spe_sub)$geo_sample_id))[1]
  
  for (n_clusters in cluster_numbers) {
    label_col <- make_label_col(n_clusters)
    labels <- colData(spe_sub)[[label_col]]
    
    for (boundary_k in boundary_k_values) {
      message("  K = ", n_clusters, "; BD k = ", boundary_k)
      
      local_col <- make_local_col(label_col, boundary_k)
      
      bd <- boundary_density(
        spatialCoords(spe_sub),
        labels,
        k = boundary_k,
        n_threads = parallel_threads,
        adjust = TRUE
      )
      
      colData(spe_out)[[local_col]] <- bd$local_boundary_density
      
      bd_summary <- rbind(
        bd_summary,
        data.frame(
          sample_id = sample_id,
          geo_sample_id = geo_sample_id,
          method = "smoothclust_sensitivity_boundary_k",
          label_source = "smoothclust",
          label_col = label_col,
          local_boundary_density_col = local_col,
          n_clusters = n_clusters,
          n_labels = length(unique(labels)),
          n_spots = ncol(spe_sub),
          boundary_k = boundary_k,
          is_default = boundary_k == default_boundary_k,
          smoothing_method = smoothing_method,
          smoothing_bandwidth = smoothing_bandwidth,
          mean_discordant = bd$mean_discordant,
          boundary_density = bd$boundary_density,
          expected_boundary_density = bd$expected_boundary_density,
          relative_boundary_density = bd$relative_boundary_density,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  spe_list[[sample_id]] <- spe_out
}

spe_plot <- do.call(cbind, spe_list)


# save outputs

parameters <- data.frame(
  dataset_id = dataset_id,
  bd_input_file = bd_input_file,
  spe_input_file = spe_input_file,
  input_boundary_k = paste(sort(unique(bd_input$boundary_k)), collapse = ","),
  clustering_scope = "per_sample",
  label_source = "smoothclust",
  smoothing_method = smoothing_method,
  smoothing_bandwidth = smoothing_bandwidth,
  boundary_k_values = paste(boundary_k_values, collapse = ","),
  default_boundary_k = default_boundary_k,
  cluster_numbers = paste(cluster_numbers, collapse = ","),
  parallel_threads = parallel_threads,
  stringsAsFactors = FALSE
)

write.csv(
  bd_summary,
  file = file.path(results_dir, "GBM_BD_sensitivity_boundary_k_summary.csv"),
  row.names = FALSE
)

write.csv(
  parameters,
  file = file.path(results_dir, "GBM_BD_sensitivity_boundary_k_parameters.csv"),
  row.names = FALSE
)

saveRDS(
  spe_plot,
  file = spe_plot_file
)
