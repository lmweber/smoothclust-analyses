################################################################################
# Script to evaluate boundary density metric
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from root directory of smoothclust-analyses repo

library(smoothclust)
library(Matrix)
library(SpatialExperiment)
library(scater)
library(scran)
library(scuttle)

dataset_id <- "gbm_greenwald_2024"

repo_dir <- "."
data_dir <- file.path(repo_dir, "data", dataset_id)
processed_dir <- file.path(data_dir, "processed")
results_dir <- file.path(repo_dir, "results", dataset_id, "boundary_density")

input_file <- file.path(processed_dir, "spe_gbm_greenwald_2024.rds")

smoothing_method <- "uniform"
smoothing_bandwidth <- 0.05
boundary_k <- 6
min_gene_total_counts <- 100
hvg_prop <- 0.1
n_pcs <- 50
cluster_numbers <- 3:8
seed <- 123


# load data

spe <- readRDS(input_file)

sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))


# multi-sample gene filtering: keep genes passing per-sample gene count
# threshold in at least one sample

gene_filter_by_sample <- sapply(sample_ids, function(sample_id) {
  ix <- as.character(colData(spe)$sample_id) == sample_id
  Matrix::rowSums(assay(spe, "counts")[, ix, drop = FALSE]) >= min_gene_total_counts
})
gene_filter <- rowSums(gene_filter_by_sample) > 0
table(gene_filter)

if ("is_mito" %in% names(rowData(spe))) {
  gene_filter <- gene_filter & !rowData(spe)$is_mito
}
table(gene_filter)

spe <- spe[gene_filter, ]

dim(spe)


# run smoothclust separately for each sample

spe_list <- lapply(sample_ids, function(sample_id) {
  message("Running smoothclust: ", sample_id)
  spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
  smoothclust(
    spe_sub,
    assay_name = "counts",
    method = smoothing_method,
    bandwidth = smoothing_bandwidth, 
    n_threads = 14
  )
})

names(spe_list) <- sample_ids
spe <- do.call(cbind, spe_list)


# normalize, select HVGs, and compute PCA

spe <- logNormCounts(spe, assay.type = "counts_smooth")

dec <- modelGeneVar(
  spe,
  block = factor(colData(spe)$sample_id)
)

top_hvgs <- getTopHVGs(dec, prop = hvg_prop)

set.seed(seed)
spe <- runPCA(
  spe,
  subset_row = top_hvgs,
  ncomponents = n_pcs
)


# cluster domains

pca <- reducedDim(spe, "PCA")
label_cols <- paste0("smoothclust_k", cluster_numbers)

for (i in seq_along(cluster_numbers)) {
  n_clusters <- cluster_numbers[i]
  label_col <- label_cols[i]
  
  set.seed(seed + n_clusters)
  clus <- kmeans(pca, centers = n_clusters, nstart = 25)$cluster
  colData(spe)[[label_col]] <- factor(clus)
}


# calculate boundary density

bd_summary <- data.frame()

for (label_col in label_cols) {
  local_col <- paste0("local_bd_", label_col)
  colData(spe)[[local_col]] <- NA_real_
  n_clusters <- as.integer(sub("smoothclust_k", "", label_col))
  
  for (sample_id in sample_ids) {
    ix <- which(as.character(colData(spe)$sample_id) == sample_id)
    
    bd <- boundary_density(
      spatialCoords(spe)[ix, , drop = FALSE],
      colData(spe)[[label_col]][ix],
      k = boundary_k,
      adjust = "analytic"
    )
    
    colData(spe)[[local_col]][ix] <- bd$local_boundary_density
    
    bd_summary <- rbind(
      bd_summary,
      data.frame(
        sample_id = sample_id,
        geo_sample_id = unique(as.character(colData(spe)$geo_sample_id[ix]))[1],
        label_col = label_col,
        n_clusters = n_clusters,
        n_spots = length(ix),
        boundary_k = boundary_k,
        smoothing_method = smoothing_method,
        smoothing_bandwidth = smoothing_bandwidth,
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


# save outputs

parameters <- data.frame(
  dataset_id = dataset_id,
  input_file = input_file,
  smoothing_method = smoothing_method,
  smoothing_bandwidth = smoothing_bandwidth,
  boundary_k = boundary_k,
  min_gene_total_counts = min_gene_total_counts,
  hvg_prop = hvg_prop,
  n_pcs = n_pcs,
  cluster_numbers = paste(cluster_numbers, collapse = ","),
  seed = seed,
  stringsAsFactors = FALSE
)

write.csv(
  bd_summary,
  file = file.path(results_dir, "gbm_boundary_density_summary.csv"),
  row.names = FALSE
)

write.csv(
  parameters,
  file = file.path(results_dir, "gbm_boundary_density_parameters.csv"),
  row.names = FALSE
)
