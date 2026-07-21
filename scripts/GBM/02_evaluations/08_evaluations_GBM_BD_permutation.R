################################################################################
# Script to evaluate permutation-adjusted boundary density
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(smoothclust)
library(BiocNeighbors)
library(SpatialExperiment)

dataset_id <- "GBM"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)

smoothclust_spe_plot_file <- file.path(results_dir, "spe_GBM_BD_plot.rds")

label_source <- "smoothclust"
cluster_numbers <- 6L
boundary_k <- 6L
n_permutations <- 1000L
seed <- 123
n_threads <- 14L

representative_sample_id <- "zh1019t1"
representative_n_clusters <- 6L

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# helper functions

permutation_boundary_density_distribution <- function(
    spatial_coords,
    labels,
    k,
    n_permutations,
    seed = NULL,
    n_threads = 1) {
  
  nn_data <- BiocNeighbors::findKNN(
    spatial_coords,
    k = k,
    get.index = TRUE,
    get.distance = FALSE,
    num.threads = n_threads
  )
  neigh <- nn_data$index
  n_neighbors <- rep.int(ncol(neigh), nrow(neigh))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  permutation_boundary_density <- numeric(n_permutations)
  for (i in seq_len(n_permutations)) {
    permuted_labels <- sample(labels)
    perm_neigh_labels <- matrix(permuted_labels[neigh], ncol = ncol(neigh))
    permutation_boundary_density[i] <- 
      mean(rowSums(permuted_labels != perm_neigh_labels) / n_neighbors)
  }
  
  permutation_boundary_density
}


# load data

spe <- readRDS(smoothclust_spe_plot_file)

sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))
label_cols <- paste0("smoothclust_k", cluster_numbers)
stopifnot(all(label_cols %in% names(colData(spe))))


# calculate permutation-adjusted BD summaries

bd_summary <- data.frame()
representative_null <- NULL

for (sample_i in seq_along(sample_ids)) {
  sample_id <- sample_ids[sample_i]
  message("Processing sample: ", sample_id)
  
  spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
  geo_sample_id <- unique(as.character(colData(spe_sub)$geo_sample_id))[1]
  
  for (n_clusters in cluster_numbers) {
    label_col <- paste0("smoothclust_k", n_clusters)
    labels <- colData(spe_sub)[[label_col]]
    permutation_seed <- seed + sample_i * 1000L + n_clusters
    
    bd <- boundary_density(
      spatialCoords(spe_sub),
      labels,
      k = boundary_k,
      n_threads = n_threads,
      adjust = "permutation",
      n_permutations = n_permutations,
      seed = permutation_seed
    )
    
    bd_summary <- rbind(
      bd_summary,
      data.frame(
        sample_id = sample_id,
        geo_sample_id = geo_sample_id,
        label_source = label_source,
        label_col = label_col,
        n_clusters = n_clusters,
        n_spots = ncol(spe_sub),
        boundary_k = boundary_k,
        n_permutations = bd$n_permutations,
        permutation_seed = permutation_seed,
        mean_discordant = bd$mean_discordant,
        boundary_density = bd$boundary_density,
        permutation_mean = bd$permutation_mean,
        permutation_sd = bd$permutation_sd,
        permutation_relative_boundary_density = 
          bd$permutation_relative_boundary_density,
        permutation_z = bd$permutation_z,
        stringsAsFactors = FALSE
      )
    )
    
    if (sample_id == representative_sample_id &&
        n_clusters == representative_n_clusters) {
      permutation_values <- permutation_boundary_density_distribution(
        spatial_coords = spatialCoords(spe_sub),
        labels = labels,
        k = boundary_k,
        n_permutations = n_permutations,
        seed = permutation_seed,
        n_threads = n_threads
      )
      
      representative_null <- data.frame(
        sample_id = sample_id,
        geo_sample_id = geo_sample_id,
        label_source = label_source,
        label_col = label_col,
        n_clusters = n_clusters,
        n_spots = ncol(spe_sub),
        boundary_k = boundary_k,
        n_permutations = n_permutations,
        permutation_seed = permutation_seed,
        permutation_id = seq_len(n_permutations),
        permutation_boundary_density = permutation_values,
        stringsAsFactors = FALSE
      )
    }
  }
}

stopifnot(!is.null(representative_null))


# save outputs

parameters <- data.frame(
  dataset_id = dataset_id,
  smoothclust_spe_plot_file = smoothclust_spe_plot_file,
  label_source = label_source,
  cluster_numbers = paste(cluster_numbers, collapse = ","),
  representative_sample_id = representative_sample_id,
  representative_n_clusters = representative_n_clusters,
  boundary_k = boundary_k,
  adjustment = "permutation",
  n_permutations = n_permutations,
  seed = seed,
  n_threads = n_threads,
  stringsAsFactors = FALSE
)

write.csv(
  bd_summary,
  file = file.path(results_dir, "GBM_BD_permutation_summary.csv"),
  row.names = FALSE
)

write.csv(
  representative_null,
  file = file.path(
    results_dir,
    paste0(
      "GBM_BD_permutation_null_",
      representative_sample_id,
      "_k",
      representative_n_clusters,
      ".csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  parameters,
  file = file.path(results_dir, "GBM_BD_permutation_parameters.csv"),
  row.names = FALSE
)
