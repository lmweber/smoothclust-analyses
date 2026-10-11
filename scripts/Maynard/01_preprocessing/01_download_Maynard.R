################################################################################
# Download the Maynard et al. human DLPFC Visium dataset
################################################################################

# Run from the smoothclust-analyses repository root directory.

dataset_id <- "Maynard"
raw_dir <- file.path("data", dataset_id, "raw")
processed_dir <- file.path("data", dataset_id, "processed")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

url <- paste0(
  "https://www.dropbox.com/s/f4wcvtdq428y73p/",
  "Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata?dl=1"
)
file <- file.path(
  raw_dir, "Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata"
)

if (!file.exists(file)) {
  part <- paste0(file, ".part")
  options(timeout = 7200)
  stopifnot(download.file(url, part, mode = "wb", method = "libcurl") == 0L)
  stopifnot(file.rename(part, file))
}

manifest <- data.frame(
  type = "author_processed_spot_level_sce",
  source = "spatialLIBD archive for Maynard et al. (2021)",
  url = url,
  file = file,
  bytes = file.info(file)$size,
  md5 = unname(tools::md5sum(file)),
  stringsAsFactors = FALSE
)
write.csv(
  manifest, file.path(processed_dir, "download_manifest.csv"),
  row.names = FALSE
)
