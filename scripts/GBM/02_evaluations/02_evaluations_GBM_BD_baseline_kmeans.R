################################################################################
# Script to evaluate non-spatial k-means baseline boundary density
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

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

baseline_method <- "baseline_unsmoothed_kmeans"
assay_name <- "counts"
normalized_assay_name <- "logcounts"
size_factor_col <- "nonmito_umi"
size_factor_reference <-
  "all_nonmitochondrial_raw_genes_before_low_expression_filtering"
size_factor_centering <- "per_sample_arithmetic_mean_one"
normalization_method <- "library_size"
normalization_pseudo_count <- 1
clustering_method <- "kmeans"
boundary_k <- 6
min_gene_total_counts <- 100
hvg_prop <- 0.1
n_pcs <- 50
cluster_numbers <- 3:8
kmeans_nstart <- 25
seed <- 123

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


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
  Matrix::rowSums(assay(spe, assay_name)[, ix, drop = FALSE]) >= min_gene_total_counts
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


# run non-spatial preprocessing and clustering separately for each sample

label_cols <- paste0("baseline_kmeans_k", cluster_numbers)
spe_list <- list()
bd_summary <- data.frame()

for (sample_id in sample_ids) {
  message("Processing sample: ", sample_id)
  spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
  size_factors <- as.numeric(colData(spe_sub)[[size_factor_col]])
  size_factors <- size_factors / mean(size_factors)
  names(size_factors) <- colnames(spe_sub)
  
  spe_sub <- logNormCounts(
    spe_sub,
    size.factors = size_factors,
    center.size.factors = FALSE,
    pseudo.count = normalization_pseudo_count,
    assay.type = assay_name,
    name = normalized_assay_name
  )
  
  dec <- modelGeneVar(spe_sub, assay.type = normalized_assay_name)
  top_hvgs <- getTopHVGs(dec, prop = hvg_prop)
  
  set.seed(seed)
  spe_sub <- runPCA(
    spe_sub,
    subset_row = top_hvgs,
    ncomponents = n_pcs,
    exprs_values = normalized_assay_name
  )
  
  pca <- reducedDim(spe_sub, "PCA")
  
  for (i in seq_along(cluster_numbers)) {
    n_clusters <- cluster_numbers[i]
    label_col <- label_cols[i]
    local_col <- paste0("local_bd_", label_col)
    
    set.seed(seed + n_clusters)
    clus <- kmeans(pca, centers = n_clusters, nstart = kmeans_nstart)$cluster
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
        method = baseline_method,
        label_col = label_col,
        n_clusters = n_clusters,
        n_spots = ncol(spe_sub),
        boundary_k = boundary_k,
        assay_name = normalized_assay_name,
        smoothing_method = NA_character_,
        smoothing_bandwidth = NA_real_,
        clustering_method = clustering_method,
        mean_discordant = bd$mean_discordant,
        boundary_density = bd$boundary_density,
        expected_boundary_density = bd$expected_boundary_density,
        relative_boundary_density = bd$relative_boundary_density,
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
  baseline_method = baseline_method,
  clustering_scope = "per_sample",
  assay_name = assay_name,
  normalized_assay_name = normalized_assay_name,
  normalization_method = normalization_method,
  normalization_pseudo_count = normalization_pseudo_count,
  size_factor_col = size_factor_col,
  size_factor_reference = size_factor_reference,
  size_factor_centering = size_factor_centering,
  clustering_method = clustering_method,
  boundary_k = boundary_k,
  min_gene_total_counts = min_gene_total_counts,
  hvg_prop = hvg_prop,
  n_pcs = n_pcs,
  cluster_numbers = paste(cluster_numbers, collapse = ","),
  kmeans_nstart = kmeans_nstart,
  seed = seed,
  stringsAsFactors = FALSE
)

write.csv(
  bd_summary,
  file = file.path(results_dir, "GBM_BD_baseline_kmeans_summary.csv"),
  row.names = FALSE
)

write.csv(
  parameters,
  file = file.path(results_dir, "GBM_BD_baseline_kmeans_parameters.csv"),
  row.names = FALSE
)


# create and save lightweight SPE object for plotting

spe_plot <- spe[0, ]
reducedDims(spe_plot) <- S4Vectors::SimpleList()

saveRDS(
  spe_plot,
  file = file.path(results_dir, "spe_GBM_baseline_kmeans_BD_plot.rds")
)
