################################################################################
# Script to evaluate boundary density for external spatial clustering methods
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

dataset_id <- "GBM"

repo_dir <- "."
data_dir <- file.path(repo_dir, "data", dataset_id)
processed_dir <- file.path(data_dir, "processed")
results_dir <- file.path(repo_dir, "results", dataset_id)

input_file <- file.path(processed_dir, "spe_GBM.rds")
spe_plot_file <- file.path(results_dir, "spe_GBM_external_clustering_BD_plot.rds")

methods_to_run <- c("bayesspace", "banksy", "clustsignal")
assay_name <- "counts"
normalized_assay_name <- "logcounts"
size_factor_col <- "nonmito_umi"
size_factor_reference <-
  "all_nonmitochondrial_raw_genes_before_low_expression_filtering"
size_factor_centering <- "per_sample_arithmetic_mean_one"
normalization_method <- "library_size"
normalization_pseudo_count <- 1
boundary_k <- 6
min_gene_total_counts <- 100
hvg_prop <- 0.1
n_pcs <- 50
cluster_numbers <- 5L
seed <- 123
parallel_threads <- 14

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


# runtime helpers

script_start_time <- Sys.time()

elapsed_seconds <- function(start_time, end_time = Sys.time()) {
  as.numeric(difftime(end_time, start_time, units = "secs"))
}

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0L) {
    NA_character_
  } else {
    paste(x, collapse = ",")
  }
}


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

store_boundary_density <- function(
    spe_out,
    labels,
    method,
    label_col,
    n_clusters_requested,
    sample_id,
    boundary_k) {
  
  local_col <- paste0("local_bd_", label_col)
  labels <- factor(labels)
  
  if (length(labels) != ncol(spe_out)) {
    stop("Label length does not match number of spots for ", label_col)
  }
  if (anyNA(labels)) {
    stop("Missing labels found for ", label_col)
  }
  
  bd <- boundary_density(
    spatialCoords(spe_out),
    labels,
    k = boundary_k,
    n_threads = parallel_threads,
    adjust = "analytic"
  )
  
  colData(spe_out)[[label_col]] <- labels
  colData(spe_out)[[local_col]] <- bd$local_boundary_density
  
  summary <- data.frame(
    sample_id = sample_id,
    geo_sample_id = unique(as.character(colData(spe_out)$geo_sample_id))[1],
    method = method,
    label_col = label_col,
    local_boundary_density_col = local_col,
    n_clusters_requested = n_clusters_requested,
    n_labels = length(unique(labels)),
    n_spots = ncol(spe_out),
    boundary_k = boundary_k,
    mean_discordant = bd$mean_discordant,
    boundary_density = bd$boundary_density,
    expected_boundary_density = bd$expected_boundary_density,
    relative_boundary_density = bd$relative_boundary_density,
    stringsAsFactors = FALSE
  )
  
  list(spe = spe_out, summary = summary)
}

run_bayesspace <- function(
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

run_banksy <- function(
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

run_clustsignal <- function(spe_sub, sample_seed, size_factors) {
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
input_n_samples <- length(sample_ids)
input_n_spots <- ncol(spe)
input_n_genes <- nrow(spe)


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
filtered_n_samples <- length(sample_ids)
filtered_n_spots <- ncol(spe)
filtered_n_genes <- nrow(spe)

dim(spe)


# run external clustering methods separately for each sample

spe_list <- list()
bd_summary <- data.frame()
runtime_by_sample <- data.frame()
analysis_start_time <- Sys.time()

for (sample_i in seq_along(sample_ids)) {
  sample_id <- sample_ids[sample_i]
  message("Processing sample: ", sample_id)
  
  sample_seed <- seed + sample_i * 1000L
  spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
  size_factors <- center_size_factors(spe_sub)
  spe_out <- spe_sub[0, ]
  
  if ("bayesspace" %in% methods_to_run) {
    message("  Running BayesSpace")
    method_start_time <- Sys.time()
    bayesspace_labels <- run_bayesspace(
      spe_sub,
      cluster_numbers,
      sample_seed,
      size_factors
    )
    for (label_col in names(bayesspace_labels)) {
      n_clusters <- as.integer(sub("bayesspace_k", "", label_col))
      res <- store_boundary_density(
        spe_out = spe_out,
        labels = bayesspace_labels[[label_col]],
        method = "bayesspace",
        label_col = label_col,
        n_clusters_requested = n_clusters,
        sample_id = sample_id,
        boundary_k = boundary_k
      )
      spe_out <- res$spe
      bd_summary <- rbind(bd_summary, res$summary)
    }
    runtime_by_sample <- rbind(
      runtime_by_sample,
      data.frame(
        dataset_id = dataset_id,
        analysis_name = "external_clustering_method_comparison",
        sample_id = sample_id,
        geo_sample_id = unique(as.character(colData(spe_sub)$geo_sample_id))[1],
        method = "bayesspace",
        method_display = "BayesSpace",
        n_clusters_requested = paste(cluster_numbers, collapse = ","),
        n_label_sets = length(bayesspace_labels),
        n_labels_returned = paste(
          vapply(
            bayesspace_labels,
            function(x) length(unique(x)),
            integer(1)
          ),
          collapse = ","
        ),
        n_spots = ncol(spe_sub),
        n_genes = nrow(spe_sub),
        boundary_k = boundary_k,
        elapsed_seconds = elapsed_seconds(method_start_time),
        requested_threads = parallel_threads,
        sample_seed = sample_seed,
        stringsAsFactors = FALSE
      )
    )
  }
  
  if ("banksy" %in% methods_to_run) {
    message("  Running Banksy")
    method_start_time <- Sys.time()
    banksy_labels <- run_banksy(
      spe_sub,
      cluster_numbers,
      sample_seed,
      size_factors
    )
    for (label_col in names(banksy_labels)) {
      n_clusters <- as.integer(sub("banksy_k", "", label_col))
      res <- store_boundary_density(
        spe_out = spe_out,
        labels = banksy_labels[[label_col]],
        method = "banksy",
        label_col = label_col,
        n_clusters_requested = n_clusters,
        sample_id = sample_id,
        boundary_k = boundary_k
      )
      spe_out <- res$spe
      bd_summary <- rbind(bd_summary, res$summary)
    }
    runtime_by_sample <- rbind(
      runtime_by_sample,
      data.frame(
        dataset_id = dataset_id,
        analysis_name = "external_clustering_method_comparison",
        sample_id = sample_id,
        geo_sample_id = unique(as.character(colData(spe_sub)$geo_sample_id))[1],
        method = "banksy",
        method_display = "Banksy",
        n_clusters_requested = paste(cluster_numbers, collapse = ","),
        n_label_sets = length(banksy_labels),
        n_labels_returned = paste(
          vapply(
            banksy_labels,
            function(x) length(unique(x)),
            integer(1)
          ),
          collapse = ","
        ),
        n_spots = ncol(spe_sub),
        n_genes = nrow(spe_sub),
        boundary_k = boundary_k,
        elapsed_seconds = elapsed_seconds(method_start_time),
        requested_threads = parallel_threads,
        sample_seed = sample_seed,
        stringsAsFactors = FALSE
      )
    )
  }
  
  if ("clustsignal" %in% methods_to_run) {
    message("  Running clustSIGNAL")
    method_start_time <- Sys.time()
    clustsignal_labels <- run_clustsignal(
      spe_sub,
      sample_seed,
      size_factors
    )
    for (label_col in names(clustsignal_labels)) {
      n_clusters <- if (grepl("^clustsignal_k[0-9]+$", label_col)) {
        as.integer(sub("clustsignal_k", "", label_col))
      } else {
        NA_integer_
      }
      res <- store_boundary_density(
        spe_out = spe_out,
        labels = clustsignal_labels[[label_col]],
        method = "clustsignal",
        label_col = label_col,
        n_clusters_requested = n_clusters,
        sample_id = sample_id,
        boundary_k = boundary_k
      )
      spe_out <- res$spe
      bd_summary <- rbind(bd_summary, res$summary)
    }
    runtime_by_sample <- rbind(
      runtime_by_sample,
      data.frame(
        dataset_id = dataset_id,
        analysis_name = "external_clustering_method_comparison",
        sample_id = sample_id,
        geo_sample_id = unique(as.character(colData(spe_sub)$geo_sample_id))[1],
        method = "clustsignal",
        method_display = "ClustSIGNAL",
        n_clusters_requested = NA_character_,
        n_label_sets = length(clustsignal_labels),
        n_labels_returned = paste(
          vapply(
            clustsignal_labels,
            function(x) length(unique(x)),
            integer(1)
          ),
          collapse = ","
        ),
        n_spots = ncol(spe_sub),
        n_genes = nrow(spe_sub),
        boundary_k = boundary_k,
        elapsed_seconds = elapsed_seconds(method_start_time),
        requested_threads = clustsignal_threads,
        sample_seed = sample_seed,
        stringsAsFactors = FALSE
      )
    )
  }
  
  spe_list[[sample_id]] <- spe_out
}

spe_plot <- do.call(cbind, spe_list)
analysis_end_time <- Sys.time()


# save outputs

parameters <- data.frame(
  dataset_id = dataset_id,
  input_file = input_file,
  methods_to_run = paste(methods_to_run, collapse = ","),
  clustering_scope = "per_sample",
  assay_name = assay_name,
  normalized_assay_name = normalized_assay_name,
  normalization_method = normalization_method,
  normalization_pseudo_count = normalization_pseudo_count,
  size_factor_col = size_factor_col,
  size_factor_reference = size_factor_reference,
  size_factor_centering = size_factor_centering,
  boundary_k = boundary_k,
  min_gene_total_counts = min_gene_total_counts,
  hvg_prop = hvg_prop,
  n_pcs = n_pcs,
  cluster_numbers = paste(cluster_numbers, collapse = ","),
  seed = seed,
  parallel_threads = parallel_threads,
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
  clustsignal_threads = clustsignal_threads,
  clustsignal_parameter_mode = clustsignal_parameter_mode,
  bayesspace_version = as.character(packageVersion("BayesSpace")),
  banksy_version = as.character(packageVersion("Banksy")),
  clustsignal_version = as.character(packageVersion("clustSIGNAL")),
  stringsAsFactors = FALSE
)

write.csv(
  bd_summary,
  file = file.path(results_dir, "GBM_BD_external_clustering_summary.csv"),
  row.names = FALSE
)

write.csv(
  parameters,
  file = file.path(results_dir, "GBM_BD_external_clustering_parameters.csv"),
  row.names = FALSE
)

saveRDS(
  spe_plot,
  file = spe_plot_file
)


# save runtime outputs

script_end_time <- Sys.time()

runtime_summary <- do.call(
  rbind,
  lapply(
    split(runtime_by_sample, runtime_by_sample$method),
    function(x) {
      data.frame(
        dataset_id = dataset_id,
        analysis_name = "external_clustering_method_comparison",
        input_file = input_file,
        method = unique(x$method)[1],
        method_display = unique(x$method_display)[1],
        n_samples_input = input_n_samples,
        n_spots_input = input_n_spots,
        n_genes_input = input_n_genes,
        n_samples_after_filter = filtered_n_samples,
        n_spots_after_filter = filtered_n_spots,
        n_genes_after_filter = filtered_n_genes,
        n_samples_timed = nrow(x),
        n_clusters_requested = collapse_unique(x$n_clusters_requested),
        n_label_sets = sum(x$n_label_sets),
        n_labels_returned = collapse_unique(x$n_labels_returned),
        boundary_k = boundary_k,
        total_elapsed_seconds = sum(x$elapsed_seconds),
        median_sample_elapsed_seconds = median(x$elapsed_seconds),
        min_sample_elapsed_seconds = min(x$elapsed_seconds),
        max_sample_elapsed_seconds = max(x$elapsed_seconds),
        requested_threads = collapse_unique(as.character(x$requested_threads)),
        methods_to_run = paste(methods_to_run, collapse = ","),
        analysis_start_time = as.character(analysis_start_time),
        analysis_end_time = as.character(analysis_end_time),
        analysis_elapsed_seconds = elapsed_seconds(
          analysis_start_time,
          analysis_end_time
        ),
        script_start_time = as.character(script_start_time),
        script_end_time = as.character(script_end_time),
        script_elapsed_seconds = elapsed_seconds(
          script_start_time,
          script_end_time
        ),
        bayesspace_version = as.character(packageVersion("BayesSpace")),
        banksy_version = as.character(packageVersion("Banksy")),
        clustsignal_version = as.character(packageVersion("clustSIGNAL")),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  runtime_by_sample,
  file = file.path(results_dir, "GBM_runtime_external_clustering_by_sample.csv"),
  row.names = FALSE
)

write.csv(
  runtime_summary,
  file = file.path(results_dir, "GBM_runtime_external_clustering_summary.csv"),
  row.names = FALSE
)
