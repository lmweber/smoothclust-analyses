################################################################################
# Script to preprocess dataset
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(VisiumIO)
library(BiocIO)
library(Matrix)
library(SpatialExperiment)

dataset_id <- "GBM"
geo_accession <- "GSE237183"
min_spot_umi <- 1000
max_spot_pct_mito <- 20

repo_dir <- "."
data_dir <- file.path(repo_dir, "data", dataset_id)
extracted_dir <- file.path(data_dir, "extracted")
processed_dir <- file.path(data_dir, "processed")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)


# find samples

matrix_files <- list.files(
  extracted_dir,
  pattern = "_matrix\\.mtx\\.gz$",
  full.names = TRUE
)

sample_prefixes <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_files))
sample_ids <- sub("^GSM[0-9]+_", "", sample_prefixes)
geo_sample_ids <- sub("_.*$", "", sample_prefixes)

sample_metadata <- data.frame(
  geo_sample_id = geo_sample_ids,
  sample_id = sample_ids,
  sample_prefix = sample_prefixes,
  stringsAsFactors = FALSE
)

sample_metadata <- sample_metadata[order(sample_metadata$sample_prefix), ]


# helper functions

flat_file <- function(sample_prefix, suffix) {
  normalizePath(file.path(extracted_dir, paste0(sample_prefix, "_", suffix)))
}

prepare_visium_folder <- function(sample_prefix, sample_dir) {
  outs_dir <- file.path(sample_dir, "outs")
  matrix_dir <- file.path(outs_dir, "filtered_feature_bc_matrix")
  spatial_dir <- file.path(outs_dir, "spatial")
  dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)
  
  matrix_suffixes <- c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
  for (suffix in matrix_suffixes) {
    from <- flat_file(sample_prefix, suffix)
    to <- file.path(matrix_dir, suffix)
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    ok <- file.symlink(from, to)
    if (!isTRUE(ok)) {
      file.copy(from, to, overwrite = TRUE)
    }
  }
  
  R.utils::gunzip(
    flat_file(sample_prefix, "tissue_positions_list.csv.gz"),
    destname = file.path(spatial_dir, "tissue_positions_list.csv"),
    remove = FALSE,
    overwrite = TRUE
  )
  R.utils::gunzip(
    flat_file(sample_prefix, "scalefactors_json.json.gz"),
    destname = file.path(spatial_dir, "scalefactors_json.json"),
    remove = FALSE,
    overwrite = TRUE
  )
  R.utils::gunzip(
    flat_file(sample_prefix, "tissue_lowres_image.png.gz"),
    destname = file.path(spatial_dir, "tissue_lowres_image.png"),
    remove = FALSE,
    overwrite = TRUE
  )
}


# load samples

tmp_dir <- tempfile("visiumio_")
dir.create(tmp_dir)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

sample_dirs <- file.path(tmp_dir, sample_metadata$sample_prefix)
mapply(prepare_visium_folder, sample_metadata$sample_prefix, sample_dirs)

vis <- TENxVisiumList(
  sampleFolders = sample_dirs,
  sample_ids = sample_metadata$sample_id,
  processing = "filtered",
  images = "lowres",
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres")
)

spe <- BiocIO::import(vis)


# filter spots and genes

spe <- spe[, as.integer(colData(spe)$in_tissue) == 1]

# Reproduce the Greenwald et al. spot-level QC: retain in-tissue spots with at
# least 1,000 total UMIs and at most 20% mitochondrial UMIs.
counts <- assay(spe, "counts")
gene_symbols <- as.character(rowData(spe)$Symbol)
is_mito <- !is.na(gene_symbols) &
  grepl("^MT-", gene_symbols, ignore.case = TRUE)
count_values <- if (inherits(counts, "sparseMatrix")) {
  counts@x
} else {
  as.vector(counts)
}

if (any(!is.finite(count_values)) || any(count_values < 0)) {
  stop("The raw count assay must contain finite, non-negative values.")
}
if (!any(!is_mito)) {
  stop("No non-mitochondrial genes were identified.")
}

sum_umi <- Matrix::colSums(counts)
mito_umi <- Matrix::colSums(counts[is_mito, , drop = FALSE])
nonmito_umi <- Matrix::colSums(counts[!is_mito, , drop = FALSE])
pct_mito <- 100 * mito_umi / sum_umi
spot_filter <- is.finite(sum_umi) &
  is.finite(pct_mito) &
  sum_umi >= min_spot_umi &
  pct_mito <= max_spot_pct_mito &
  nonmito_umi > 0

qc_sample_id <- as.character(colData(spe)$sample_id)
qc_summary <- do.call(rbind, lapply(sort(unique(qc_sample_id)), function(sample_id) {
  ix <- qc_sample_id == sample_id
  data.frame(
    sample_id = sample_id,
    n_in_tissue = sum(ix),
    n_retained = sum(ix & spot_filter),
    n_removed = sum(ix & !spot_filter),
    n_below_min_umi = sum(ix & sum_umi < min_spot_umi),
    n_above_max_pct_mito = sum(ix & pct_mito > max_spot_pct_mito, na.rm = TRUE),
    n_nonfinite_or_zero_umi = sum(ix & (!is.finite(pct_mito) | sum_umi <= 0)),
    stringsAsFactors = FALSE
  )
}))

spe <- spe[, spot_filter]

# Remove genes with zero counts across the retained spots.
spe <- spe[Matrix::rowSums(assay(spe, "counts")) > 0, ]


# add QC metrics on the retained object

counts <- assay(spe, "counts")
gene_symbols <- as.character(rowData(spe)$Symbol)
is_mito <- !is.na(gene_symbols) &
  grepl("^MT-", gene_symbols, ignore.case = TRUE)

metadata_match <- match(as.character(colData(spe)$sample_id), sample_metadata$sample_id)
colData(spe)$barcode <- colnames(spe)
colData(spe)$dataset_id <- dataset_id
colData(spe)$geo_accession <- geo_accession
colData(spe)$geo_sample_id <- sample_metadata$geo_sample_id[metadata_match]
colData(spe)$sample_prefix <- sample_metadata$sample_prefix[metadata_match]
colnames(spe) <- paste(colData(spe)$sample_id, colData(spe)$barcode, sep = "_")

rowData(spe)$symbol <- gene_symbols
rowData(spe)$is_mito <- is_mito
rowData(spe)$total_counts <- Matrix::rowSums(counts)
rowData(spe)$detected_spots <- Matrix::rowSums(counts > 0)

colData(spe)$sum_umi <- Matrix::colSums(counts)
colData(spe)$detected_genes <- Matrix::colSums(counts > 0)
colData(spe)$mito_umi <- Matrix::colSums(counts[is_mito, , drop = FALSE])
colData(spe)$nonmito_umi <- Matrix::colSums(counts[!is_mito, , drop = FALSE])
colData(spe)$pct_mito <- 100 * colData(spe)$mito_umi / colData(spe)$sum_umi

stopifnot(
  all(colData(spe)$sum_umi >= min_spot_umi),
  all(colData(spe)$pct_mito <= max_spot_pct_mito),
  all(is.finite(colData(spe)$nonmito_umi)),
  all(colData(spe)$nonmito_umi > 0),
  all(colData(spe)$sum_umi ==
        colData(spe)$mito_umi + colData(spe)$nonmito_umi)
)


# save outputs

metadata(spe)$dataset_id <- dataset_id
metadata(spe)$geo_accession <- geo_accession
metadata(spe)$sample_metadata <- sample_metadata
metadata(spe)$spot_qc <- list(
  min_spot_umi = min_spot_umi,
  max_spot_pct_mito = max_spot_pct_mito,
  summary = qc_summary
)

saveRDS(
  spe,
  file = file.path(processed_dir, "spe_GBM.rds")
)

write.csv(
  qc_summary,
  file = file.path(processed_dir, "GBM_preprocessing_QC_summary.csv"),
  row.names = FALSE
)
