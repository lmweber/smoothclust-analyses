################################################################################
# Script to compare spatial clustering methods against author-derived labels
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(smoothclust)
library(Matrix)
library(SpatialExperiment)
library(scater)
library(scran)
library(scuttle)
library(BayesSpace)
library(Banksy)
library(clustSIGNAL)
library(mclust)

dataset_id <- "GBM"

repo_dir <- "."
data_dir <- file.path(repo_dir, "data", dataset_id)
processed_dir <- file.path(data_dir, "processed")
results_dir <- file.path(repo_dir, "results", dataset_id)

input_file <- file.path(processed_dir, "spe_GBM.rds")
author_metadata_file <- file.path(data_dir, "raw", "Inputs", "general", "visium_metadata.csv")
spe_plot_file <- file.path(results_dir, "spe_GBM_author_label_clustering_comparison_plot.rds")

methods_to_run <- c("smoothclust", "bayesspace", "banksy", "clustsignal")
fixed_cluster_numbers <- c(3L, 5L)
author_target_cols <- c("org1", "org2")
assay_name <- "counts"
normalized_assay_name <- "logcounts"
smooth_assay_name <- "logcounts_smooth"
size_factor_col <- "nonmito_umi"
size_factor_reference <-
  "all_nonmitochondrial_raw_genes_before_low_expression_filtering"
size_factor_centering <- "per_sample_arithmetic_mean_one"
normalization_method <- "library_size"
normalization_pseudo_count <- 1
smoothclust_normalize <- TRUE
smoothclust_log <- TRUE
boundary_k <- 6
min_gene_total_counts <- 100
hvg_prop <- 0.1
n_pcs <- 50
seed <- 123
parallel_threads <- 14

# smoothclust settings
smoothclust_smoothing_method <- "uniform"
smoothclust_smoothing_bandwidth <- 0.05
smoothclust_kmeans_nstart <- 25

# BayesSpace settings
bayesspace_platform <- "Visium"
bayesspace_n_pcs <- 15
bayesspace_nrep <- 50000
bayesspace_burn_in <- 1000
bayesspace_thin <- 100
bayesspace_init_method <- "kmeans"
bayesspace_model <- "t"
bayesspace_precision <- "equal"

# Banksy settings
banksy_assay_name <- normalized_assay_name
banksy_m <- 0
banksy_lambda <- 0.8
banksy_k_geom <- 15
banksy_compute_agf <- FALSE
banksy_algo <- "kmeans"

# clustSIGNAL settings
clustsignal_threads <- parallel_threads
clustsignal_parameter_mode <- "default"

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# metadata

author_target_metadata <- data.frame(
  target_label_col = c("org1", "org2"),
  target_label_display = c("org1 zones (3 classes)", "org2 zones (5 classes)"),
  matched_n_clusters = c(3L, 5L),
  stringsAsFactors = FALSE
)

author_sample_map <- data.frame(
  sample_id = c(
    "bwh23", "bwh24", "bwh25", "bwh28", "bwh35",
    "mgh258", "mgh259",
    "zh1007inf", "zh1007nec", "zh1019inf", "zh1019t1",
    "zh8811a", "zh8811b", "zh8812", "zh881inf", "zh881t1",
    "zh916bulk", "zh916inf", "zh916t1"
  ),
  author_sample_id = c(
    "bwh23", "bwh24", "bwh25", "bwh28", "bwh35",
    "MGH258", "mgh259",
    "ZH1007inf", "ZH1007nec", "ZH1019inf", "ZH1019T1",
    "ZH8811Abulk", "ZH8811Bbulk", "ZH8812bulk", "ZH881inf", "ZH881T1",
    "ZH916bulk", "ZH916inf", "ZH916T1"
  ),
  stringsAsFactors = FALSE
)

author_sample_groups <- data.frame(
  sample_id = c(
    "mgh258",
    "zh1007inf", "zh1007nec",
    "zh1019inf", "zh1019t1",
    "zh8811a", "zh8811b", "zh8812", "zh881inf", "zh881t1",
    "zh916bulk", "zh916inf", "zh916t1"
  ),
  tumor_group = c(
    "MGH258",
    "ZH1007", "ZH1007",
    "ZH1019", "ZH1019",
    "ZH881", "ZH881", "ZH881", "ZH881", "ZH881",
    "ZH916", "ZH916", "ZH916"
  ),
  region_group = c(
    "Tumor, unqualified",
    "Infiltrating", "Necrotic",
    "Infiltrating", "T1 contrast-enhancing",
    "Bulk tumor", "Bulk tumor", "Bulk tumor", "Infiltrating", "T1 contrast-enhancing",
    "Bulk tumor", "Infiltrating", "T1 contrast-enhancing"
  ),
  stringsAsFactors = FALSE
)

analysis_sample_ids <- author_sample_groups$sample_id


# helper functions

select_hvgs <- function(spe, prop, assay_type = normalized_assay_name) {
  dec <- modelGeneVar(spe, assay.type = assay_type)
  getTopHVGs(dec, prop = prop)
}

parallel_bpparam <- function(n_threads) {
  if (.Platform$OS.type == "unix") {
    BiocParallel::MulticoreParam(workers = n_threads)
  } else {
    BiocParallel::SnowParam(workers = n_threads)
  }
}

safe_n_pcs <- function(spe, requested_n_pcs) {
  max_pcs <- min(nrow(spe), ncol(spe)) - 1
  max(1L, min(requested_n_pcs, max_pcs))
}

center_size_factors <- function(spe) {
  size_factors <- as.numeric(colData(spe)[[size_factor_col]])
  if (length(size_factors) != ncol(spe) || any(!is.finite(size_factors)) ||
      any(size_factors <= 0)) {
    stop("Invalid size factors in colData(spe)$", size_factor_col, ".")
  }
  size_factors <- size_factors / mean(size_factors)
  names(size_factors) <- colnames(spe)
  size_factors
}

is_matched_granularity <- function(method, n_clusters_requested, target_label_col) {
  if (method == "clustsignal") {
    return(TRUE)
  }
  target_ix <- match(target_label_col, author_target_metadata$target_label_col)
  !is.na(n_clusters_requested) &&
    n_clusters_requested == author_target_metadata$matched_n_clusters[target_ix]
}

attach_author_labels <- function(spe, author_metadata, author_sample_map, author_cols) {
  sample_match <- match(as.character(colData(spe)$sample_id), author_sample_map$sample_id)
  colData(spe)$author_sample_id <- author_sample_map$author_sample_id[sample_match]
  
  author_key <- paste(author_metadata$sample, author_metadata$spot_id, sep = "___")
  stopifnot(!anyDuplicated(author_key))
  
  spe_key <- paste(colData(spe)$author_sample_id, colData(spe)$barcode, sep = "___")
  author_match <- match(spe_key, author_key)
  colData(spe)$author_metadata_match <- !is.na(author_match)
  
  for (author_col in author_cols) {
    colData(spe)[[paste0("author_", author_col)]] <- author_metadata[[author_col]][author_match]
  }
  
  spe
}

store_method_results <- function(
    spe_out,
    labels,
    method,
    method_display,
    label_col,
    n_clusters_requested,
    sample_id,
    sample_seed,
    boundary_k) {
  
  labels <- factor(labels)
  if (length(labels) != ncol(spe_out)) {
    stop("Label length does not match number of spots for ", label_col)
  }
  if (anyNA(labels)) {
    stop("Missing labels found for ", label_col)
  }
  
  local_col <- paste0("local_bd_", label_col)
  bd <- boundary_density(
    spatialCoords(spe_out),
    labels,
    k = boundary_k,
    n_threads = parallel_threads,
    adjust = "analytic"
  )
  
  colData(spe_out)[[label_col]] <- labels
  colData(spe_out)[[local_col]] <- bd$local_boundary_density
  
  group_ix <- match(sample_id, author_sample_groups$sample_id)
  summary <- data.frame()
  
  for (target_i in seq_len(nrow(author_target_metadata))) {
    target_label_col <- author_target_metadata$target_label_col[target_i]
    target_col <- paste0("author_", target_label_col)
    target <- colData(spe_out)[[target_col]]
    compare_ix <- !is.na(target)
    
    if (sum(compare_ix) == 0L) {
      ari <- NA_real_
      n_target_labels <- NA_integer_
    } else {
      ari <- mclust::adjustedRandIndex(labels[compare_ix], factor(target[compare_ix]))
      n_target_labels <- length(unique(target[compare_ix]))
    }
    
    summary <- rbind(
      summary,
      data.frame(
        sample_id = sample_id,
        geo_sample_id = unique(as.character(colData(spe_out)$geo_sample_id))[1],
        author_sample_id = unique(as.character(colData(spe_out)$author_sample_id))[1],
        tumor_group = author_sample_groups$tumor_group[group_ix],
        region_group = author_sample_groups$region_group[group_ix],
        method = method,
        method_display = method_display,
        label_col = label_col,
        local_boundary_density_col = local_col,
        n_clusters_requested = n_clusters_requested,
        n_labels = length(unique(labels)),
        target_label_col = target_label_col,
        target_label_display = author_target_metadata$target_label_display[target_i],
        target_matched_n_clusters = author_target_metadata$matched_n_clusters[target_i],
        n_target_labels = n_target_labels,
        is_matched_granularity = is_matched_granularity(
          method,
          n_clusters_requested,
          target_label_col
        ),
        n_spots = ncol(spe_out),
        n_compared_spots = sum(compare_ix),
        boundary_k = boundary_k,
        mean_discordant = bd$mean_discordant,
        boundary_density = bd$boundary_density,
        expected_boundary_density = bd$expected_boundary_density,
        relative_boundary_density = bd$relative_boundary_density,
        adjusted_rand_index = ari,
        sample_seed = sample_seed,
        stringsAsFactors = FALSE
      )
    )
  }
  
  list(spe = spe_out, summary = summary)
}

run_smoothclust_labels <- function(
    spe_sub,
    cluster_numbers,
    sample_seed,
    size_factors) {
  spe_sc <- smoothclust(
    spe_sub,
    assay_name = assay_name,
    method = smoothclust_smoothing_method,
    bandwidth = smoothclust_smoothing_bandwidth,
    n_threads = parallel_threads,
    normalize = smoothclust_normalize,
    log = smoothclust_log,
    size_factors = size_factors,
    pseudo_count = normalization_pseudo_count,
    output_assay_name = smooth_assay_name
  )
  
  top_hvgs <- select_hvgs(
    spe_sc,
    hvg_prop,
    assay_type = smooth_assay_name
  )
  n_pcs_use <- safe_n_pcs(spe_sc[top_hvgs, ], n_pcs)
  
  set.seed(sample_seed)
  spe_sc <- runPCA(
    spe_sc,
    subset_row = top_hvgs,
    ncomponents = n_pcs_use,
    exprs_values = smooth_assay_name
  )
  pca <- reducedDim(spe_sc, "PCA")
  
  labels <- list()
  for (n_clusters in cluster_numbers) {
    set.seed(sample_seed + n_clusters)
    clus <- kmeans(
      pca,
      centers = n_clusters,
      nstart = smoothclust_kmeans_nstart
    )$cluster
    labels[[paste0("smoothclust_k", n_clusters)]] <- factor(clus, levels = seq_len(n_clusters))
  }
  
  labels
}

run_bayesspace_labels <- function(
    spe_sub,
    cluster_numbers,
    sample_seed,
    size_factors) {
  n_hvgs <- max(1L, min(nrow(spe_sub), ceiling(hvg_prop * nrow(spe_sub))))
  n_pcs_use <- safe_n_pcs(spe_sub, bayesspace_n_pcs)
  
  spe_bs <- spe_sub
  SingleCellExperiment::sizeFactors(spe_bs) <- size_factors
  colData(spe_bs)$row <- colData(spe_bs)$array_row
  colData(spe_bs)$col <- colData(spe_bs)$array_col
  
  spe_bs <- BayesSpace::spatialPreprocess(
    spe_bs,
    platform = bayesspace_platform,
    n.PCs = n_pcs_use,
    n.HVGs = n_hvgs,
    log.normalize = TRUE,
    assay.type = normalized_assay_name,
    BPPARAM = parallel_bpparam(parallel_threads)
  )
  
  labels <- list()
  for (n_clusters in cluster_numbers) {
    set.seed(sample_seed + n_clusters)
    spe_clustered <- BayesSpace::spatialCluster(
      spe_bs,
      q = n_clusters,
      use.dimred = "PCA",
      d = n_pcs_use,
      platform = bayesspace_platform,
      init.method = bayesspace_init_method,
      model = bayesspace_model,
      precision = bayesspace_precision,
      nrep = bayesspace_nrep,
      burn.in = bayesspace_burn_in,
      thin = bayesspace_thin
    )
    labels[[paste0("bayesspace_k", n_clusters)]] <- 
      factor(colData(spe_clustered)$spatial.cluster)
  }
  
  labels
}

run_banksy_labels <- function(
    spe_sub,
    cluster_numbers,
    sample_seed,
    size_factors) {
  spe_banksy <- logNormCounts(
    spe_sub,
    size.factors = size_factors,
    center.size.factors = FALSE,
    pseudo.count = normalization_pseudo_count,
    assay.type = assay_name,
    name = normalized_assay_name
  )
  top_hvgs <- select_hvgs(
    spe_banksy,
    hvg_prop,
    assay_type = normalized_assay_name
  )
  spe_banksy <- spe_banksy[top_hvgs, ]
  n_pcs_use <- safe_n_pcs(spe_banksy, n_pcs)
  
  spe_banksy <- Banksy::computeBanksy(
    spe_banksy,
    assay_name = banksy_assay_name,
    compute_agf = banksy_compute_agf,
    k_geom = banksy_k_geom,
    M = banksy_m,
    seed = sample_seed,
    parallel = parallel_threads > 1L,
    num_cores = parallel_threads,
    verbose = FALSE
  )
  
  spe_banksy <- Banksy::runBanksyPCA(
    spe_banksy,
    M = banksy_m,
    lambda = banksy_lambda,
    npcs = n_pcs_use,
    assay_name = banksy_assay_name,
    seed = sample_seed
  )
  
  spe_banksy <- Banksy::clusterBanksy(
    spe_banksy,
    M = banksy_m,
    lambda = banksy_lambda,
    use_pcs = TRUE,
    npcs = n_pcs_use,
    algo = banksy_algo,
    kmeans.centers = cluster_numbers,
    seed = sample_seed
  )
  
  labels <- list()
  banksy_cols <- Banksy::clusterNames(spe_banksy)
  for (n_clusters in cluster_numbers) {
    col_pattern <- paste0("kmeans", n_clusters, "$")
    banksy_col <- banksy_cols[grepl(col_pattern, banksy_cols)]
    if (length(banksy_col) != 1L) {
      stop("Could not identify Banksy cluster column for K = ", n_clusters)
    }
    labels[[paste0("banksy_k", n_clusters)]] <- factor(colData(spe_banksy)[[banksy_col]])
  }
  
  labels
}

run_clustsignal_labels <- function(spe_sub, sample_seed, size_factors) {
  spe_cs <- logNormCounts(
    spe_sub,
    size.factors = size_factors,
    center.size.factors = FALSE,
    pseudo.count = normalization_pseudo_count,
    assay.type = assay_name,
    name = normalized_assay_name
  )
  top_hvgs <- select_hvgs(
    spe_cs,
    hvg_prop,
    assay_type = normalized_assay_name
  )
  spe_cs <- spe_cs[top_hvgs, ]
  
  set.seed(sample_seed)
  res <- clustSIGNAL::clustSIGNAL(
    spe_cs,
    samples = "sample_id",
    threads = clustsignal_threads,
    outputs = "s"
  )
  
  list(clustsignal_default = factor(colData(res$spe_final)$ClustSIGNAL))
}


# load and filter data

spe <- readRDS(input_file)
author_metadata <- read.csv(author_metadata_file, stringsAsFactors = FALSE)

stopifnot(all(c("sample", "spot_id", author_target_cols) %in% names(author_metadata)))
stopifnot(all(analysis_sample_ids %in% unique(as.character(colData(spe)$sample_id))))
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

spe <- spe[, as.character(colData(spe)$sample_id) %in% analysis_sample_ids]
spe <- attach_author_labels(spe, author_metadata, author_sample_map, author_target_cols)

sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))

# multi-sample gene filtering among the annotated samples only
gene_filter_by_sample <- sapply(sample_ids, function(sample_id) {
  ix <- as.character(colData(spe)$sample_id) == sample_id
  Matrix::rowSums(assay(spe, assay_name)[, ix, drop = FALSE]) >=
    min_gene_total_counts
})
gene_filter <- rowSums(gene_filter_by_sample) > 0

gene_filter <- gene_filter & !rowData(spe)$is_mito

spe <- spe[gene_filter, ]

# remove spots with zero counts after gene filtering
spot_counts_after_gene_filter <- Matrix::colSums(assay(spe, assay_name))
spot_filter <- spot_counts_after_gene_filter > 0
spe <- spe[, spot_filter]
sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))


# run methods separately for each sample

spe_list <- list()
comparison_summary <- data.frame()

for (sample_i in seq_along(sample_ids)) {
  sample_id <- sample_ids[sample_i]
  message("Processing sample: ", sample_id)
  
  sample_seed <- seed + sample_i * 1000L
  spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
  size_factors <- center_size_factors(spe_sub)
  spe_out <- spe_sub[0, ]
  
  if ("smoothclust" %in% methods_to_run) {
    message("  Running smoothclust")
    smoothclust_labels <- run_smoothclust_labels(
      spe_sub,
      fixed_cluster_numbers,
      sample_seed,
      size_factors
    )
    for (label_col in names(smoothclust_labels)) {
      n_clusters <- as.integer(sub("smoothclust_k", "", label_col))
      res <- store_method_results(
        spe_out = spe_out,
        labels = smoothclust_labels[[label_col]],
        method = "smoothclust",
        method_display = "smoothclust",
        label_col = label_col,
        n_clusters_requested = n_clusters,
        sample_id = sample_id,
        sample_seed = sample_seed,
        boundary_k = boundary_k
      )
      spe_out <- res$spe
      comparison_summary <- rbind(comparison_summary, res$summary)
    }
  }
  
  if ("bayesspace" %in% methods_to_run) {
    message("  Running BayesSpace")
    bayesspace_labels <- run_bayesspace_labels(
      spe_sub,
      fixed_cluster_numbers,
      sample_seed,
      size_factors
    )
    for (label_col in names(bayesspace_labels)) {
      n_clusters <- as.integer(sub("bayesspace_k", "", label_col))
      res <- store_method_results(
        spe_out = spe_out,
        labels = bayesspace_labels[[label_col]],
        method = "bayesspace",
        method_display = "BayesSpace",
        label_col = label_col,
        n_clusters_requested = n_clusters,
        sample_id = sample_id,
        sample_seed = sample_seed,
        boundary_k = boundary_k
      )
      spe_out <- res$spe
      comparison_summary <- rbind(comparison_summary, res$summary)
    }
  }
  
  if ("banksy" %in% methods_to_run) {
    message("  Running Banksy")
    banksy_labels <- run_banksy_labels(
      spe_sub,
      fixed_cluster_numbers,
      sample_seed,
      size_factors
    )
    for (label_col in names(banksy_labels)) {
      n_clusters <- as.integer(sub("banksy_k", "", label_col))
      res <- store_method_results(
        spe_out = spe_out,
        labels = banksy_labels[[label_col]],
        method = "banksy",
        method_display = "Banksy",
        label_col = label_col,
        n_clusters_requested = n_clusters,
        sample_id = sample_id,
        sample_seed = sample_seed,
        boundary_k = boundary_k
      )
      spe_out <- res$spe
      comparison_summary <- rbind(comparison_summary, res$summary)
    }
  }
  
  if ("clustsignal" %in% methods_to_run) {
    message("  Running clustSIGNAL")
    clustsignal_labels <- run_clustsignal_labels(
      spe_sub,
      sample_seed,
      size_factors
    )
    for (label_col in names(clustsignal_labels)) {
      res <- store_method_results(
        spe_out = spe_out,
        labels = clustsignal_labels[[label_col]],
        method = "clustsignal",
        method_display = "ClustSIGNAL",
        label_col = label_col,
        n_clusters_requested = NA_integer_,
        sample_id = sample_id,
        sample_seed = sample_seed,
        boundary_k = boundary_k
      )
      spe_out <- res$spe
      comparison_summary <- rbind(comparison_summary, res$summary)
    }
  }
  
  spe_list[[sample_id]] <- spe_out
}

spe_plot <- do.call(cbind, spe_list)


# save outputs

parameters <- data.frame(
  dataset_id = dataset_id,
  input_file = input_file,
  author_metadata_file = author_metadata_file,
  methods_to_run = paste(methods_to_run, collapse = ","),
  analysis_sample_ids = paste(analysis_sample_ids, collapse = ","),
  author_target_cols = paste(author_target_cols, collapse = ","),
  clustering_scope = "per_author_labeled_sample",
  assay_name = assay_name,
  normalized_assay_name = normalized_assay_name,
  smooth_assay_name = smooth_assay_name,
  normalization_method = normalization_method,
  normalization_pseudo_count = normalization_pseudo_count,
  smoothclust_normalize = smoothclust_normalize,
  smoothclust_log = smoothclust_log,
  size_factor_col = size_factor_col,
  size_factor_reference = size_factor_reference,
  size_factor_centering = size_factor_centering,
  boundary_k = boundary_k,
  min_gene_total_counts = min_gene_total_counts,
  hvg_prop = hvg_prop,
  n_pcs = n_pcs,
  fixed_cluster_numbers = paste(fixed_cluster_numbers, collapse = ","),
  seed = seed,
  parallel_threads = parallel_threads,
  smoothclust_smoothing_method = smoothclust_smoothing_method,
  smoothclust_smoothing_bandwidth = smoothclust_smoothing_bandwidth,
  bayesspace_platform = bayesspace_platform,
  bayesspace_n_pcs = bayesspace_n_pcs,
  bayesspace_nrep = bayesspace_nrep,
  bayesspace_burn_in = bayesspace_burn_in,
  bayesspace_thin = bayesspace_thin,
  bayesspace_init_method = bayesspace_init_method,
  bayesspace_model = bayesspace_model,
  bayesspace_precision = bayesspace_precision,
  banksy_assay_name = banksy_assay_name,
  banksy_m = banksy_m,
  banksy_lambda = banksy_lambda,
  banksy_k_geom = banksy_k_geom,
  banksy_compute_agf = banksy_compute_agf,
  banksy_algo = banksy_algo,
  clustsignal_parameter_mode = clustsignal_parameter_mode,
  smoothclust_version = as.character(packageVersion("smoothclust")),
  bayesspace_version = as.character(packageVersion("BayesSpace")),
  banksy_version = as.character(packageVersion("Banksy")),
  clustsignal_version = as.character(packageVersion("clustSIGNAL")),
  mclust_version = as.character(packageVersion("mclust")),
  stringsAsFactors = FALSE
)

write.csv(
  comparison_summary,
  file = file.path(results_dir, "GBM_BD_author_label_clustering_comparison_summary.csv"),
  row.names = FALSE
)

write.csv(
  parameters,
  file = file.path(results_dir, "GBM_BD_author_label_clustering_comparison_parameters.csv"),
  row.names = FALSE
)

saveRDS(
  spe_plot,
  file = spe_plot_file
)
