################################################################################
# Preprocess the Maynard et al. human DLPFC Visium dataset
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(spatialLIBD)

dataset_id <- "Maynard"
raw_file <- file.path(
  "data", dataset_id, "raw",
  "Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata"
)
processed_dir <- file.path("data", dataset_id, "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

load(raw_file)
sce <- spatialLIBD:::.update_sce(sce)
meta <- as.data.frame(colData(sce))

spots <- data.frame(
  spot_id = meta$key,
  barcode = meta$barcode,
  dataset_id = dataset_id,
  section_id = as.character(meta$sample_name),
  sample_id = as.character(meta$sample_name),
  patient_id = meta$subject,
  position_um = meta$position,
  adjacent_replicate = meta$replicate,
  pair_id = meta$subject_position,
  array_row = meta$row,
  array_col = meta$col,
  image_row = meta$imagerow,
  image_col = meta$imagecol,
  x = meta$col / 2,
  y = meta$row * sqrt(3) / 2,
  author_layer = as.character(meta$layer_guess_reordered_short),
  author_discard = meta$discard,
  stringsAsFactors = FALSE
)

sample_metadata <- unique(spots[c(
  "section_id", "sample_id", "patient_id", "position_um",
  "adjacent_replicate", "pair_id"
)])
sample_metadata <- sample_metadata[order(
  sample_metadata$patient_id, sample_metadata$position_um,
  sample_metadata$adjacent_replicate
), ]
sample_metadata$spatialLIBD_version <-
  as.character(packageVersion("spatialLIBD"))
sample_metadata$n_spots <- as.integer(table(spots$section_id)[
  sample_metadata$section_id
])
sample_metadata$n_labeled_spots <- as.integer(tapply(
  !is.na(spots$author_layer), spots$section_id, sum
)[sample_metadata$section_id])
sample_metadata$labeled_fraction <-
  sample_metadata$n_labeled_spots / sample_metadata$n_spots
sample_metadata$n_author_discard <- as.integer(tapply(
  spots$author_discard, spots$section_id, sum
)[sample_metadata$section_id])

stopifnot(
  nrow(spots) == 47681L,
  nrow(sample_metadata) == 12L,
  sum(!is.na(spots$author_layer)) == 47329L
)

saveRDS(spots, file.path(processed_dir, "Maynard_spot_metadata.rds"))
write.csv(
  sample_metadata, file.path(processed_dir, "Maynard_sample_metadata.csv"),
  row.names = FALSE
)
