################################################################################
# Preprocess the Miyoshi et al. human Alzheimer disease Visium dataset
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(Matrix)
library(SpatialExperiment)

dataset_id <- "Miyoshi"
raw_dir <- file.path("data", dataset_id, "raw")
processed_dir <- file.path("data", dataset_id, "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# GEO compressed an already gzip-compressed RDS; two gzip connections are needed.
outer <- gzfile(file.path(
  raw_dir, "GSE233208_Human_visium_ADDS_seurat_processed.rds.gz"
), "rb")
inner <- gzcon(outer)
object <- readRDS(inner)
close(inner)

# Extract deposited Seurat slots without loading Seurat.
meta <- attr(object, "meta.data")

assays <- attr(object, "assays")
counts <- attr(assays[["Spatial"]], "counts")
images <- attr(object, "images")
coord_list <- lapply(names(images), function(image_id) {
  x <- attr(images[[image_id]], "coordinates")
  data.frame(
    spot_id = rownames(x), author_image_id = image_id,
    array_row = x$row, array_col = x$col,
    image_row = x$imagerow, image_col = x$imagecol,
    stringsAsFactors = FALSE
  )
})
coordinates <- do.call(rbind, coord_list)
# Align image coordinates with the metadata and count-matrix column order.
coordinates <- coordinates[match(rownames(meta), coordinates$spot_id), ]
stopifnot(
  identical(rownames(meta), colnames(counts)),
  identical(rownames(meta), coordinates$spot_id)
)

cortical_region <- as.character(meta$annotation)
spots <- data.frame(
  spot_id = rownames(meta),
  barcode = sub("_[0-9]+_[0-9]+$", "", rownames(meta)),
  dataset_id = dataset_id,
  section_id = as.character(meta$Sample),
  sample_id = as.character(meta$Sample),
  author_image_id = coordinates$author_image_id,
  array_row = coordinates$array_row,
  array_col = coordinates$array_col,
  image_row = coordinates$image_row,
  image_col = coordinates$image_col,
  author_cortical_region = cortical_region,
  author_network_region = ifelse(grepl("^WM", cortical_region), "WM", cortical_region),
  stringsAsFactors = FALSE,
  row.names = rownames(meta)
)

sample_data <- unique(spots[c(
  "section_id", "sample_id", "author_image_id"
)])
sample_data$n_spots <- as.integer(table(spots$section_id)[
  sample_data$section_id
])
stopifnot(
  nrow(sample_data) == 39L,
  sum(sample_data$n_spots) == 115451L,
  !anyNA(spots$author_cortical_region)
)

# Transform Visium array positions to an equilateral lattice for neighbor distances.
spatial_coords <- cbind(
  x = spots$array_col / 2,
  y = spots$array_row * sqrt(3) / 2
)
rownames(spatial_coords) <- spots$spot_id
gene_data <- DataFrame(
  gene_name = rownames(counts),
  symbol = rownames(counts),
  is_mito = grepl("^MT-", rownames(counts)),
  row.names = rownames(counts)
)
spe <- SpatialExperiment(
  assays = list(counts = counts),
  rowData = gene_data,
  colData = DataFrame(spots),
  spatialCoords = spatial_coords
)

metadata(spe)$dataset_id <- dataset_id
metadata(spe)$geo_accession <- "GSE233208"
metadata(spe)$source <- "Miyoshi et al. (2024), doi:10.1038/s41588-024-01961-x"
metadata(spe)$annotation_description <- paste(
  "author_cortical_region is the deposited nine-level, marker- and",
  "location-annotated BayesSpace q=9 assignment; author_network_region",
  "collapses WM1-WM3 to WM as in the authors' seven regional networks"
)
metadata(spe)$coordinate_description <- paste(
  "deposited Visium array positions transformed to an equilateral hexagonal",
  "lattice; deposited pixel positions are retained in colData"
)
metadata(spe)$sample_metadata <- sample_data

rm(object, assays, images, meta, counts)
gc(FALSE)
saveRDS(spe, file.path(processed_dir, "spe_Miyoshi.rds"))
write.csv(
  sample_data, file.path(processed_dir, "Miyoshi_sample_metadata.csv"),
  row.names = FALSE
)
