################################################################################
# Script to evaluate smoothing method sensitivity for boundary density
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(smoothclust)
library(Matrix)
library(SpatialExperiment)
library(scater)
library(scran)
library(mclust)

dataset_id <- "GBM"

repo_dir <- "."
data_dir <- file.path(repo_dir, "data", dataset_id)
processed_dir <- file.path(data_dir, "processed")
results_dir <- file.path(repo_dir, "results", dataset_id)

input_file <- file.path(processed_dir, "spe_GBM.rds")
spe_plot_file <- file.path(
  results_dir,
  "spe_GBM_BD_sensitivity_smoothing_method_plot.rds"
)

smoothing_methods <- c("uniform", "kernel", "knn")
assay_name <- "counts"
smooth_assay_name <- "logcounts_smooth"
size_factor_col <- "nonmito_umi"
size_factor_reference <-
  "all_nonmitochondrial_raw_genes_before_low_expression_filtering"
size_factor_centering <- "per_sample_arithmetic_mean_one"
smoothclust_normalize <- TRUE
smoothclust_log <- TRUE
smoothclust_pseudo_count <- 1
smoothing_bandwidth <- 0.05
smoothing_k <- 18
boundary_k <- 6
min_gene_total_counts <- 100
hvg_prop <- 0.1
n_pcs <- 50
cluster_numbers <- 3:8
kmeans_nstart <- 25
seed <- 123
parallel_threads <- 14

method_pairs <- data.frame(
  method_x = c("uniform", "uniform", "kernel"),
  method_y = c("kernel", "knn", "knn"),
  stringsAsFactors = FALSE
)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# helper functions

smooth_by_method <- function(spe_sub, method, size_factors) {
  if (method %in% c("uniform", "kernel")) {
    smoothclust(
      spe_sub,
      assay_name = assay_name,
      method = method,
      bandwidth = smoothing_bandwidth,
      n_threads = parallel_threads,
      normalize = smoothclust_normalize,
      log = smoothclust_log,
      size_factors = size_factors,
      pseudo_count = smoothclust_pseudo_count,
      output_assay_name = smooth_assay_name
    )
  } else if (method == "knn") {
    smoothclust(
      spe_sub,
      assay_name = assay_name,
      method = method,
      k = smoothing_k,
      n_threads = parallel_threads,
      normalize = smoothclust_normalize,
      log = smoothclust_log,
      size_factors = size_factors,
      pseudo_count = smoothclust_pseudo_count,
      output_assay_name = smooth_assay_name
    )
  } else {
    stop("Unknown smoothing method: ", method)
  }
}

make_label_col <- function(method, n_clusters) {
  paste0("smoothclust_", method, "_k", n_clusters)
}


# load data

spe <- readRDS(input_file)

if (!size_factor_col %in% names(colData(spe))) {
  stop("Missing required column in colData(spe): ", size_factor_col)
}
nonmito_umi <- colData(spe)[[size_factor_col]]
if (!is.numeric(nonmito_umi) || length(nonmito_umi) != ncol(spe) ||
    any(!is.finite(nonmito_umi)) || any(nonmito_umi <= 0)) {
  stop("colData(spe)$", size_factor_col,
       " must contain one positive finite numeric value per spot.")
}
if (!"is_mito" %in% names(rowData(spe)) ||
    !is.logical(rowData(spe)$is_mito) ||
    length(rowData(spe)$is_mito) != nrow(spe) ||
    anyNA(rowData(spe)$is_mito)) {
  stop("rowData(spe)$is_mito must contain one non-missing logical value per gene.")
}

sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))


# multi-sample gene filtering: keep genes passing per-sample gene count
# threshold in at least one sample

gene_filter_by_sample <- sapply(sample_ids, function(sample_id) {
  ix <- as.character(colData(spe)$sample_id) == sample_id
  Matrix::rowSums(assay(spe, assay_name)[, ix, drop = FALSE]) >= 
    min_gene_total_counts
})
gene_filter <- rowSums(gene_filter_by_sample) > 0
table(gene_filter)

gene_filter <- gene_filter & !rowData(spe)$is_mito
table(gene_filter)

spe <- spe[gene_filter, ]


# remove spots with zero counts after gene filtering

spot_counts_after_gene_filter <- Matrix::colSums(assay(spe, assay_name))
spot_filter <- spot_counts_after_gene_filter > 0
spe <- spe[, spot_filter]
sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))

dim(spe)


# run smoothing methods and clustering separately for each sample

spe_list <- list()
bd_summary <- data.frame()
ari_summary <- data.frame()

for (sample_id in sample_ids) {
  message("Processing sample: ", sample_id)
  
  spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
  size_factors <- as.numeric(colData(spe_sub)[[size_factor_col]])
  size_factors <- size_factors / mean(size_factors)
  names(size_factors) <- colnames(spe_sub)
  spe_out <- spe_sub[0, ]
  geo_sample_id <- unique(as.character(colData(spe_sub)$geo_sample_id))[1]
  labels_by_method <- list()
  
  for (method_i in seq_along(smoothing_methods)) {
    smoothing_method <- smoothing_methods[method_i]
    message("  Smoothing method: ", smoothing_method)
    
    spe_method <- smooth_by_method(
      spe_sub,
      smoothing_method,
      size_factors
    )
    
    dec <- modelGeneVar(spe_method, assay.type = smooth_assay_name)
    top_hvgs <- getTopHVGs(dec, prop = hvg_prop)
    n_pcs_use <- min(n_pcs, length(top_hvgs), ncol(spe_method) - 1L)
    
    # Match the primary GBM analysis and reset the seed for each method so
    # stochastic steps are controlled across smoothing-method comparisons.
    set.seed(seed)
    spe_method <- runPCA(
      spe_method,
      subset_row = top_hvgs,
      ncomponents = n_pcs_use,
      exprs_values = smooth_assay_name
    )
    
    pca <- reducedDim(spe_method, "PCA")
    labels_by_method[[smoothing_method]] <- list()
    
    for (n_clusters in cluster_numbers) {
      label_col <- make_label_col(smoothing_method, n_clusters)
      local_col <- paste0("local_bd_", label_col)
      
      set.seed(seed + n_clusters)
      clus <- kmeans(
        pca,
        centers = n_clusters,
        nstart = kmeans_nstart
      )$cluster
      labels <- factor(clus, levels = seq_len(n_clusters))
      
      bd <- boundary_density(
        spatialCoords(spe_method),
        labels,
        k = boundary_k,
        n_threads = parallel_threads,
        adjust = "analytic"
      )
      
      colData(spe_out)[[label_col]] <- labels
      colData(spe_out)[[local_col]] <- bd$local_boundary_density
      labels_by_method[[smoothing_method]][[as.character(n_clusters)]] <- labels
      
      bd_summary <- rbind(
        bd_summary,
        data.frame(
          sample_id = sample_id,
          geo_sample_id = geo_sample_id,
          method = "smoothclust_sensitivity_smoothing_method",
          label_col = label_col,
          local_boundary_density_col = local_col,
          n_clusters = n_clusters,
          n_labels = length(unique(labels)),
          n_spots = ncol(spe_method),
          boundary_k = boundary_k,
          assay_name = smooth_assay_name,
          smoothing_method = smoothing_method,
          smoothing_bandwidth = if (smoothing_method %in% c("uniform", "kernel")) {
            smoothing_bandwidth
          } else {
            NA_real_
          },
          smoothing_k = if (smoothing_method == "knn") smoothing_k else NA_integer_,
          hvg_prop = hvg_prop,
          n_hvgs = length(top_hvgs),
          n_pcs = n_pcs_use,
          clustering_method = "kmeans",
          kmeans_nstart = kmeans_nstart,
          mean_discordant = bd$mean_discordant,
          boundary_density = bd$boundary_density,
          expected_boundary_density = bd$expected_boundary_density,
          relative_boundary_density = bd$relative_boundary_density,
          excess_boundary_density = bd$excess_boundary_density,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  for (n_clusters in cluster_numbers) {
    for (pair_i in seq_len(nrow(method_pairs))) {
      method_x <- method_pairs$method_x[pair_i]
      method_y <- method_pairs$method_y[pair_i]
      label_col_x <- make_label_col(method_x, n_clusters)
      label_col_y <- make_label_col(method_y, n_clusters)
      labels_x <- labels_by_method[[method_x]][[as.character(n_clusters)]]
      labels_y <- labels_by_method[[method_y]][[as.character(n_clusters)]]
      
      ari_summary <- rbind(
        ari_summary,
        data.frame(
          sample_id = sample_id,
          geo_sample_id = geo_sample_id,
          method_x = method_x,
          method_y = method_y,
          method_pair = paste(method_x, method_y, sep = "_vs_"),
          label_col_x = label_col_x,
          label_col_y = label_col_y,
          n_clusters = n_clusters,
          n_spots = ncol(spe_sub),
          ari = adjustedRandIndex(labels_x, labels_y),
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
  input_file = input_file,
  smoothing_methods = paste(smoothing_methods, collapse = ","),
  method_pairs = paste(
    paste(method_pairs$method_x, method_pairs$method_y, sep = " vs "),
    collapse = ","
  ),
  clustering_scope = "per_sample",
  assay_name = assay_name,
  smooth_assay_name = smooth_assay_name,
  size_factor_col = size_factor_col,
  size_factor_reference = size_factor_reference,
  size_factor_centering = size_factor_centering,
  smoothclust_normalize = smoothclust_normalize,
  smoothclust_log = smoothclust_log,
  smoothclust_pseudo_count = smoothclust_pseudo_count,
  smoothing_bandwidth = smoothing_bandwidth,
  smoothing_k = smoothing_k,
  boundary_k = boundary_k,
  min_gene_total_counts = min_gene_total_counts,
  hvg_prop = hvg_prop,
  n_pcs = n_pcs,
  cluster_numbers = paste(cluster_numbers, collapse = ","),
  kmeans_nstart = kmeans_nstart,
  seed = seed,
  pca_seed_strategy = "seed for each sample and smoothing method",
  kmeans_seed_strategy = "seed + n_clusters for each sample and smoothing method",
  parallel_threads = parallel_threads,
  stringsAsFactors = FALSE
)

write.csv(
  bd_summary,
  file = file.path(results_dir, "GBM_BD_sensitivity_smoothing_method_summary.csv"),
  row.names = FALSE
)

write.csv(
  ari_summary,
  file = file.path(results_dir, "GBM_BD_sensitivity_smoothing_method_ari.csv"),
  row.names = FALSE
)

write.csv(
  parameters,
  file = file.path(results_dir, "GBM_BD_sensitivity_smoothing_method_parameters.csv"),
  row.names = FALSE
)

saveRDS(
  spe_plot,
  file = spe_plot_file
)
