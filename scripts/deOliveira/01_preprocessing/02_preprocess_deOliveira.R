################################################################################
# Preprocess illustrative de Oliveira Visium HD and Xenium sections
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(arrow)
library(DropletUtils)
library(Matrix)
library(SpatialExperiment)
library(SpatialExperimentIO)

dataset_id <- "deOliveira"
raw_dir <- file.path("data", dataset_id, "raw")
processed_dir <- file.path("data", dataset_id, "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 123
sketch_n <- 25000
min_detected_fraction <- 0.001

# Rank genes by variance of log-normalized expression on a reproducible sketch.
select_variable_genes <- function(counts_mat, library_size, keep, sketch, n) {
  size_factors <- library_size[sketch] / mean(library_size[sketch])
  logcounts <- counts_mat[keep, sketch, drop = FALSE] %*%
    Diagonal(x = 1 / size_factors)
  logcounts <- log1p(logcounts) / log(2)
  gene_mean <- Matrix::rowMeans(logcounts)
  gene_var <- Matrix::rowMeans(logcounts ^ 2) - gene_mean ^ 2
  names(gene_var) <- rownames(logcounts)
  names(sort(gene_var, decreasing = TRUE))[seq_len(min(n, length(gene_var)))]
}


# Visium HD: full filtered 8-micron P2 CRC section

hd_dir <- file.path(
  raw_dir,
  "VisiumHD_HumanColon_Oliveira",
  "binned_outputs",
  "square_008um"
)
hd_sce <- read10xCounts(
  file.path(hd_dir, "filtered_feature_bc_matrix.h5"),
  col.names = TRUE
)
hd_counts <- counts(hd_sce)
# Preserve full-library sizes before selecting genes for downstream normalization.
hd_library_size <- Matrix::colSums(hd_counts)
hd_detected <- Matrix::colSums(hd_counts > 0)
hd_keep <- hd_library_size > 0
hd_sce <- hd_sce[, hd_keep]
hd_counts <- counts(hd_sce)
hd_library_size <- hd_library_size[hd_keep]
hd_detected <- hd_detected[hd_keep]

hd_positions <- read_parquet(
  file.path(hd_dir, "spatial", "tissue_positions.parquet")
)
hd_clustering <- read.csv(file.path(hd_dir, "clustering.csv.gz"))
hd_deconvolution <- read.csv(file.path(hd_dir, "deconvolution.csv.gz"))
# Match annotation rows by barcode because input tables may use different orders.
position_match <- match(colnames(hd_sce), hd_positions$barcode)
cluster_match <- match(colnames(hd_sce), hd_clustering$barcode)
decon_match <- match(colnames(hd_sce), hd_deconvolution$barcode)
stopifnot(!anyNA(position_match), !anyNA(cluster_match))

hd_symbol <- rowData(hd_sce)$Symbol
hd_gene_detected <- Matrix::rowSums(hd_counts > 0)
hd_gene_keep <- !grepl("^MT-", hd_symbol, ignore.case = TRUE) &
  hd_gene_detected >= ceiling(min_detected_fraction * ncol(hd_sce))
set.seed(seed)
hd_sketch <- sample(seq_len(ncol(hd_sce)), min(sketch_n, ncol(hd_sce)))
hd_hvgs <- select_variable_genes(
  hd_counts,
  hd_library_size,
  hd_gene_keep,
  hd_sketch,
  n = 300
)

hd_coldata <- DataFrame(
  barcode = colnames(hd_sce),
  sample_id = "P2CRC",
  platform = "Visium HD",
  bin_size_um = 8,
  array_row = hd_positions$array_row[position_match],
  array_col = hd_positions$array_col[position_match],
  pxl_row_in_fullres = hd_positions$pxl_row_in_fullres[position_match],
  pxl_col_in_fullres = hd_positions$pxl_col_in_fullres[position_match],
  library_size = hd_library_size,
  n_detected = hd_detected,
  author_cluster = factor(hd_clustering$cluster[cluster_match]),
  decon_class = hd_deconvolution$DeconClass[decon_match],
  decon_label_1 = hd_deconvolution$DeconLabel1[decon_match],
  decon_label_2 = hd_deconvolution$DeconLabel2[decon_match],
  row.names = colnames(hd_sce)
)
hd_spe <- SpatialExperiment(
  assays = list(
    counts = as(hd_counts[hd_hvgs, , drop = FALSE], "dgCMatrix")
  ),
  rowData = rowData(hd_sce)[hd_hvgs, , drop = FALSE],
  colData = hd_coldata,
  # Visium HD coordinates are square-bin indices; one unit represents 8 um.
  spatialCoords = cbind(
    x = hd_coldata$array_col,
    y = hd_coldata$array_row
  )
)
saveRDS(
  hd_spe,
  file.path(processed_dir, "spe_deOliveira_VisiumHD_P2CRC_8um.rds")
)


# Xenium: full non-empty P2 CRC cell matrix

xenium_dir <- file.path(raw_dir, "Xenium_HumanColon_Oliveira")
xenium_spe <- readXeniumSXE(xenium_dir, addParquetPaths = FALSE)
# Preserve full-library sizes before selecting genes for downstream normalization.
xenium_library_size <- Matrix::colSums(counts(xenium_spe))
xenium_detected <- Matrix::colSums(counts(xenium_spe) > 0)
xenium_keep <- xenium_library_size > 0
xenium_spe <- xenium_spe[, xenium_keep]
xenium_library_size <- xenium_library_size[xenium_keep]
xenium_detected <- xenium_detected[xenium_keep]
colData(xenium_spe)$osta_sample_id <- colData(xenium_spe)$sample_id
colData(xenium_spe)$sample_id <- "P2CRC"
colData(xenium_spe)$platform <- "Xenium"
colData(xenium_spe)$library_size <- xenium_library_size
colData(xenium_spe)$n_detected <- xenium_detected

xenium_gene_detected <- Matrix::rowSums(counts(xenium_spe) > 0)
xenium_gene_keep <-
  xenium_gene_detected >= ceiling(min_detected_fraction * ncol(xenium_spe))
set.seed(seed)
xenium_sketch <- sample(
  seq_len(ncol(xenium_spe)),
  min(sketch_n, ncol(xenium_spe))
)
xenium_hvgs <- select_variable_genes(
  counts(xenium_spe),
  xenium_library_size,
  xenium_gene_keep,
  xenium_sketch,
  n = 300
)
xenium_spe <- xenium_spe[xenium_hvgs, ]
saveRDS(
  xenium_spe,
  file.path(processed_dir, "spe_deOliveira_Xenium_P2CRC.rds")
)


# Preprocessing summary

preprocessing_summary <- rbind(
  data.frame(
    platform = "Visium HD",
    sample_id = "P2CRC",
    spatial_unit = "8-micron square bin",
    n_input_locations = length(hd_keep),
    n_retained_locations = ncol(hd_spe),
    n_input_genes = nrow(hd_sce),
    n_selected_genes = nrow(hd_spe),
    n_zero_selected_gene_counts = sum(
      Matrix::colSums(counts(hd_spe)) == 0
    ),
    median_library_size = median(hd_coldata$library_size),
    median_detected_genes = median(hd_coldata$n_detected),
    n_author_cluster_labels = sum(!is.na(hd_coldata$author_cluster)),
    n_author_deconvolution_labels = sum(!is.na(hd_coldata$decon_label_1))
  ),
  data.frame(
    platform = "Xenium",
    sample_id = "P2CRC",
    spatial_unit = "segmented cell",
    n_input_locations = length(xenium_keep),
    n_retained_locations = ncol(xenium_spe),
    n_input_genes = length(xenium_gene_detected),
    n_selected_genes = nrow(xenium_spe),
    n_zero_selected_gene_counts = sum(
      Matrix::colSums(counts(xenium_spe)) == 0
    ),
    median_library_size = median(colData(xenium_spe)$library_size),
    median_detected_genes = median(colData(xenium_spe)$n_detected),
    n_author_cluster_labels = NA_integer_,
    n_author_deconvolution_labels = NA_integer_
  )
)
write.csv(
  preprocessing_summary,
  file.path(processed_dir, "deOliveira_preprocessing_summary.csv"),
  row.names = FALSE
)
