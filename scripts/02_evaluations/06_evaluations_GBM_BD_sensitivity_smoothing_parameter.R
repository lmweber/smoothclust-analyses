################################################################################
# Script to evaluate smoothing parameter sensitivity for boundary density
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(smoothclust)
library(Matrix)
library(SpatialExperiment)
library(scater)
library(scran)

dataset_id <- "GBM"

repo_dir <- "."
data_dir <- file.path(repo_dir, "data", dataset_id)
processed_dir <- file.path(data_dir, "processed")
results_dir <- file.path(repo_dir, "results", dataset_id)

input_file <- file.path(processed_dir, "spe_GBM.rds")
spe_plot_file <- file.path(
  results_dir,
  "spe_GBM_BD_sensitivity_smoothing_parameter_plot.rds"
)

assay_name <- "counts"
smooth_assay_name <- "logcounts_smooth"
size_factor_col <- "nonmito_umi"
size_factor_reference <-
  "all_nonmitochondrial_raw_genes_before_low_expression_filtering"
size_factor_centering <- "per_sample_arithmetic_mean_one"
smoothclust_normalize <- TRUE
smoothclust_log <- TRUE
smoothclust_pseudo_count <- 1
uniform_bandwidths <- c(0.010, 0.025, 0.050, 0.075, 0.100, 0.150)
kernel_bandwidths <- c(0.010, 0.025, 0.050, 0.075, 0.100, 0.150)
knn_smoothing_k <- c(2L, 6L, 12L, 18L, 30L, 42L)
default_smoothing_bandwidth <- 0.050
default_smoothing_k <- 18L
boundary_k <- 6
min_gene_total_counts <- 100
hvg_prop <- 0.1
n_pcs <- 50
cluster_numbers <- 3:8
kmeans_nstart <- 25
seed <- 123
parallel_threads <- 14

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# helper functions

format_parameter_id <- function(parameter_type, value) {
  if (parameter_type == "bandwidth") {
    gsub("\\.", "p", sprintf("b%.3f", value))
  } else if (parameter_type == "smoothing_k") {
    sprintf("ks%03d", as.integer(value))
  } else {
    stop("Unknown smoothing parameter type: ", parameter_type)
  }
}

make_smoothing_grid <- function() {
  smoothing_grid <- rbind(
    data.frame(
      smoothing_method = "uniform",
      smoothing_parameter_type = "bandwidth",
      smoothing_parameter_value = uniform_bandwidths,
      smoothing_bandwidth = uniform_bandwidths,
      smoothing_k = NA_integer_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      smoothing_method = "kernel",
      smoothing_parameter_type = "bandwidth",
      smoothing_parameter_value = kernel_bandwidths,
      smoothing_bandwidth = kernel_bandwidths,
      smoothing_k = NA_integer_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      smoothing_method = "knn",
      smoothing_parameter_type = "smoothing_k",
      smoothing_parameter_value = knn_smoothing_k,
      smoothing_bandwidth = NA_real_,
      smoothing_k = knn_smoothing_k,
      stringsAsFactors = FALSE
    )
  )
  
  smoothing_grid$smoothing_setting_id <- paste(
    smoothing_grid$smoothing_method,
    mapply(
      format_parameter_id,
      smoothing_grid$smoothing_parameter_type,
      smoothing_grid$smoothing_parameter_value
    ),
    sep = "_"
  )
  
  smoothing_grid$is_default <- with(
    smoothing_grid,
    (
      smoothing_method %in% c("uniform", "kernel") &
        smoothing_bandwidth == default_smoothing_bandwidth
    ) |
      (
        smoothing_method == "knn" &
          smoothing_k == default_smoothing_k
      )
  )
  
  smoothing_grid
}

smooth_by_setting <- function(spe_sub, setting, size_factors) {
  if (setting$smoothing_method %in% c("uniform", "kernel")) {
    smoothclust(
      spe_sub,
      assay_name = assay_name,
      method = setting$smoothing_method,
      bandwidth = setting$smoothing_bandwidth,
      n_threads = parallel_threads,
      normalize = smoothclust_normalize,
      log = smoothclust_log,
      size_factors = size_factors,
      pseudo_count = smoothclust_pseudo_count,
      output_assay_name = smooth_assay_name
    )
  } else if (setting$smoothing_method == "knn") {
    smoothclust(
      spe_sub,
      assay_name = assay_name,
      method = setting$smoothing_method,
      k = setting$smoothing_k,
      n_threads = parallel_threads,
      normalize = smoothclust_normalize,
      log = smoothclust_log,
      size_factors = size_factors,
      pseudo_count = smoothclust_pseudo_count,
      output_assay_name = smooth_assay_name
    )
  } else {
    stop("Unknown smoothing method: ", setting$smoothing_method)
  }
}

make_label_col <- function(smoothing_setting_id, n_clusters) {
  paste0("smoothclust_", smoothing_setting_id, "_k", n_clusters)
}

calculate_rank_stability <- function(bd_summary) {
  rank_stability <- data.frame()
  
  for (smoothing_method in unique(bd_summary$smoothing_method)) {
    for (n_clusters in sort(unique(bd_summary$n_clusters))) {
      method_k <- bd_summary[
        bd_summary$smoothing_method == smoothing_method &
          bd_summary$n_clusters == n_clusters,
      ]
      default_rows <- method_k[
        method_k$is_default,
        c("sample_id", "boundary_density")
      ]
      names(default_rows)[names(default_rows) == "boundary_density"] <-
        "default_boundary_density"
      
      for (setting_id in unique(method_k$smoothing_setting_id)) {
        setting_rows <- method_k[method_k$smoothing_setting_id == setting_id, ]
        paired_rows <- merge(
          setting_rows,
          default_rows,
          by = "sample_id",
          all = FALSE
        )
        
        rank_stability <- rbind(
          rank_stability,
          data.frame(
            smoothing_method = smoothing_method,
            smoothing_setting_id = setting_id,
            smoothing_parameter_type = unique(setting_rows$smoothing_parameter_type),
            smoothing_parameter_value = unique(setting_rows$smoothing_parameter_value),
            smoothing_bandwidth = unique(setting_rows$smoothing_bandwidth),
            smoothing_k = unique(setting_rows$smoothing_k),
            is_default = unique(setting_rows$is_default),
            n_clusters = n_clusters,
            n_samples = nrow(paired_rows),
            spearman_rho = cor(
              paired_rows$boundary_density,
              paired_rows$default_boundary_density,
              method = "spearman"
            ),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  rank_stability
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


# run smoothing parameter grid and clustering separately for each sample

smoothing_grid <- make_smoothing_grid()
spe_list <- list()
bd_summary <- data.frame()

for (sample_i in seq_along(sample_ids)) {
  sample_id <- sample_ids[sample_i]
  message("Processing sample: ", sample_id)
  
  sample_seed <- seed + sample_i * 1000L
  spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
  size_factors <- as.numeric(colData(spe_sub)[[size_factor_col]])
  size_factors <- size_factors / mean(size_factors)
  names(size_factors) <- colnames(spe_sub)
  spe_out <- spe_sub[0, ]
  geo_sample_id <- unique(as.character(colData(spe_sub)$geo_sample_id))[1]
  
  for (setting_i in seq_len(nrow(smoothing_grid))) {
    setting <- smoothing_grid[setting_i, ]
    message("  Smoothing setting: ", setting$smoothing_setting_id)
    
    spe_setting <- smooth_by_setting(
      spe_sub,
      setting,
      size_factors
    )
    
    dec <- modelGeneVar(spe_setting, assay.type = smooth_assay_name)
    top_hvgs <- getTopHVGs(dec, prop = hvg_prop)
    n_pcs_use <- min(n_pcs, length(top_hvgs), ncol(spe_setting) - 1L)
    
    set.seed(sample_seed)
    spe_setting <- runPCA(
      spe_setting,
      subset_row = top_hvgs,
      ncomponents = n_pcs_use,
      exprs_values = smooth_assay_name
    )
    
    pca <- reducedDim(spe_setting, "PCA")
    
    for (n_clusters in cluster_numbers) {
      label_col <- make_label_col(setting$smoothing_setting_id, n_clusters)
      local_col <- paste0("local_bd_", label_col)
      
      set.seed(sample_seed + n_clusters)
      clus <- kmeans(
        pca,
        centers = n_clusters,
        nstart = kmeans_nstart
      )$cluster
      labels <- factor(clus, levels = seq_len(n_clusters))
      
      bd <- boundary_density(
        spatialCoords(spe_setting),
        labels,
        k = boundary_k,
        n_threads = parallel_threads,
        adjust = "analytic"
      )
      
      colData(spe_out)[[label_col]] <- labels
      colData(spe_out)[[local_col]] <- bd$local_boundary_density
      
      bd_summary <- rbind(
        bd_summary,
        data.frame(
          sample_id = sample_id,
          geo_sample_id = geo_sample_id,
          method = "smoothclust_sensitivity_smoothing_parameter",
          label_col = label_col,
          local_boundary_density_col = local_col,
          n_clusters = n_clusters,
          n_labels = length(unique(labels)),
          n_spots = ncol(spe_setting),
          boundary_k = boundary_k,
          assay_name = smooth_assay_name,
          smoothing_method = setting$smoothing_method,
          smoothing_setting_id = setting$smoothing_setting_id,
          smoothing_parameter_type = setting$smoothing_parameter_type,
          smoothing_parameter_value = setting$smoothing_parameter_value,
          smoothing_bandwidth = setting$smoothing_bandwidth,
          smoothing_k = setting$smoothing_k,
          is_default = setting$is_default,
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
  
  spe_list[[sample_id]] <- spe_out
}

spe_plot <- do.call(cbind, spe_list)
rank_stability <- calculate_rank_stability(bd_summary)


# save outputs

parameters <- data.frame(
  dataset_id = dataset_id,
  input_file = input_file,
  clustering_scope = "per_sample",
  assay_name = assay_name,
  smooth_assay_name = smooth_assay_name,
  size_factor_col = size_factor_col,
  size_factor_reference = size_factor_reference,
  size_factor_centering = size_factor_centering,
  smoothclust_normalize = smoothclust_normalize,
  smoothclust_log = smoothclust_log,
  smoothclust_pseudo_count = smoothclust_pseudo_count,
  uniform_bandwidths = paste(uniform_bandwidths, collapse = ","),
  kernel_bandwidths = paste(kernel_bandwidths, collapse = ","),
  knn_smoothing_k = paste(knn_smoothing_k, collapse = ","),
  default_smoothing_bandwidth = default_smoothing_bandwidth,
  default_smoothing_k = default_smoothing_k,
  boundary_k = boundary_k,
  min_gene_total_counts = min_gene_total_counts,
  hvg_prop = hvg_prop,
  n_pcs = n_pcs,
  cluster_numbers = paste(cluster_numbers, collapse = ","),
  kmeans_nstart = kmeans_nstart,
  seed = seed,
  parallel_threads = parallel_threads,
  stringsAsFactors = FALSE
)

write.csv(
  bd_summary,
  file = file.path(results_dir, "GBM_BD_sensitivity_smoothing_parameter_summary.csv"),
  row.names = FALSE
)

write.csv(
  rank_stability,
  file = file.path(results_dir, "GBM_BD_sensitivity_smoothing_parameter_rank_stability.csv"),
  row.names = FALSE
)

write.csv(
  parameters,
  file = file.path(results_dir, "GBM_BD_sensitivity_smoothing_parameter_parameters.csv"),
  row.names = FALSE
)

saveRDS(
  spe_plot,
  file = spe_plot_file
)
