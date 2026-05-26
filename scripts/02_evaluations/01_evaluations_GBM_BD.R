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

dataset_id <- "GBM"

repo_dir <- "."
data_dir <- file.path(repo_dir, "data", dataset_id)
processed_dir <- file.path(data_dir, "processed")
results_dir <- file.path(repo_dir, "results", dataset_id)

input_file <- file.path(processed_dir, "spe_GBM.rds")

smoothing_method <- "uniform"
smoothing_bandwidth <- 0.05
boundary_k <- 6
min_gene_total_counts <- 100
hvg_prop <- 0.1
n_pcs <- 50
cluster_numbers <- 3:8
seed <- 123

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


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


# remove spots with zero counts after gene filtering

spot_counts_after_gene_filter <- Matrix::colSums(assay(spe, "counts"))
spot_filter <- spot_counts_after_gene_filter > 0
spe <- spe[, spot_filter]
sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))

dim(spe)


# run smoothclust and clustering separately for each sample

label_cols <- paste0("smoothclust_k", cluster_numbers)
spe_list <- list()
bd_summary <- data.frame()

for (sample_id in sample_ids) {
  message("Processing sample: ", sample_id)
  spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
  
  spe_sub <- smoothclust(
    spe_sub,
    assay_name = "counts",
    method = smoothing_method,
    bandwidth = smoothing_bandwidth, 
    n_threads = 14
  )
  
  spe_sub <- logNormCounts(spe_sub, assay.type = "counts_smooth")
  
  dec <- modelGeneVar(spe_sub)
  top_hvgs <- getTopHVGs(dec, prop = hvg_prop)
  
  set.seed(seed)
  spe_sub <- runPCA(
    spe_sub,
    subset_row = top_hvgs,
    ncomponents = n_pcs
  )
  
  pca <- reducedDim(spe_sub, "PCA")

  for (i in seq_along(cluster_numbers)) {
    n_clusters <- cluster_numbers[i]
    label_col <- label_cols[i]
    local_col <- paste0("local_bd_", label_col)
    
    set.seed(seed + n_clusters)
    clus <- kmeans(pca, centers = n_clusters, nstart = 25)$cluster
    print(table(clus))
    colData(spe_sub)[[label_col]] <- factor(clus, levels = seq_len(n_clusters))
    
    bd <- boundary_density(
      spatialCoords(spe_sub),
      colData(spe_sub)[[label_col]],
      k = boundary_k,
      adjust = "analytic"
    )
    
    colData(spe_sub)[[local_col]] <- bd$local_boundary_density
    
    bd_summary <- rbind(
      bd_summary,
      data.frame(
        sample_id = sample_id,
        geo_sample_id = unique(as.character(colData(spe_sub)$geo_sample_id))[1],
        label_col = label_col,
        n_clusters = n_clusters,
        n_spots = ncol(spe_sub),
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
  
  spe_list[[sample_id]] <- spe_sub
}

spe <- do.call(cbind, spe_list)


# save outputs

parameters <- data.frame(
  dataset_id = dataset_id,
  input_file = input_file,
  clustering_scope = "per_sample",
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
  file = file.path(results_dir, "GBM_BD_summary.csv"),
  row.names = FALSE
)

write.csv(
  parameters,
  file = file.path(results_dir, "GBM_BD_parameters.csv"),
  row.names = FALSE
)


# create and save lightweight SPE object for plotting

spe_plot <- spe[0, ]
reducedDims(spe_plot) <- S4Vectors::SimpleList()

saveRDS(
  spe_plot,
  file = file.path(results_dir, "spe_GBM_BD_plot.rds")
)
