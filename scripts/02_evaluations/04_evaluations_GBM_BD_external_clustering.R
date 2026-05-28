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
banksy_assay_name <- "logcounts"
banksy_m <- 0
banksy_lambda <- 0.8
banksy_k_geom <- 15
banksy_compute_agf <- FALSE
banksy_algo <- "kmeans"

# clustSIGNAL settings
clustsignal_threads <- parallel_threads
clustsignal_parameter_mode <- "default"

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# helper functions

select_hvgs <- function(spe, prop) {
  dec <- modelGeneVar(spe)
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
    excess_boundary_density = bd$excess_boundary_density,
    stringsAsFactors = FALSE
  )
  
  list(spe = spe_out, summary = summary)
}

run_bayesspace <- function(spe_sub, cluster_numbers, sample_seed) {
  n_hvgs <- max(1L, min(nrow(spe_sub), ceiling(hvg_prop * nrow(spe_sub))))
  n_pcs_use <- safe_n_pcs(spe_sub, bayesspace_n_pcs)
  
  spe_bs <- spe_sub
  colData(spe_bs)$row <- colData(spe_bs)$array_row
  colData(spe_bs)$col <- colData(spe_bs)$array_col
  
  spe_bs <- BayesSpace::spatialPreprocess(
    spe_bs,
    platform = bayesspace_platform,
    n.PCs = n_pcs_use,
    n.HVGs = n_hvgs,
    log.normalize = TRUE,
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

run_banksy <- function(spe_sub, cluster_numbers, sample_seed) {
  spe_banksy <- logNormCounts(spe_sub, assay.type = assay_name)
  top_hvgs <- select_hvgs(spe_banksy, hvg_prop)
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

run_clustsignal <- function(spe_sub, sample_seed) {
  spe_cs <- logNormCounts(spe_sub, assay.type = assay_name)
  top_hvgs <- select_hvgs(spe_cs, hvg_prop)
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

if ("is_mito" %in% names(rowData(spe))) {
  gene_filter <- gene_filter & !rowData(spe)$is_mito
}
table(gene_filter)

spe <- spe[gene_filter, ]


# remove spots with zero counts after gene filtering

spot_counts_after_gene_filter <- Matrix::colSums(assay(spe, assay_name))
spot_filter <- spot_counts_after_gene_filter > 0
spe <- spe[, spot_filter]
sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))

dim(spe)


# run external clustering methods separately for each sample

spe_list <- list()
bd_summary <- data.frame()

for (sample_i in seq_along(sample_ids)) {
  sample_id <- sample_ids[sample_i]
  message("Processing sample: ", sample_id)
  
  sample_seed <- seed + sample_i * 1000L
  spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
  spe_out <- spe_sub[0, ]
  
  if ("bayesspace" %in% methods_to_run) {
    message("  Running BayesSpace")
    bayesspace_labels <- run_bayesspace(spe_sub, cluster_numbers, sample_seed)
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
  }
  
  if ("banksy" %in% methods_to_run) {
    message("  Running Banksy")
    banksy_labels <- run_banksy(spe_sub, cluster_numbers, sample_seed)
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
  }
  
  if ("clustsignal" %in% methods_to_run) {
    message("  Running clustSIGNAL")
    clustsignal_labels <- run_clustsignal(spe_sub, sample_seed)
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
  }
  
  spe_list[[sample_id]] <- spe_out
}

spe_plot <- do.call(cbind, spe_list)


# save outputs

parameters <- data.frame(
  dataset_id = dataset_id,
  input_file = input_file,
  methods_to_run = paste(methods_to_run, collapse = ","),
  clustering_scope = "per_sample",
  assay_name = assay_name,
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
