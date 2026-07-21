################################################################################
# Download the Schabitz et al. inflammatory skin disease Visium dataset
################################################################################

# Run from the smoothclust-analyses repository root directory.

dataset_id <- "skin"
geo_accession <- "GSE206391"

raw_dir <- file.path(".", "data", dataset_id, "raw")
processed_dir <- file.path(".", "data", dataset_id, "processed")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

geo_files <- c(
  "GSE206391_Unfiltered_data.h5",
  "GSE206391_Preprocessed_data.h5"
)
geo_base_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/",
  geo_accession,
  "/suppl/"
)
local_files <- file.path(raw_dir, geo_files)

options(timeout = 7200)
for (i in seq_along(geo_files)) {
  download.file(
    paste0(geo_base_url, geo_files[i]),
    local_files[i],
    mode = "wb",
    method = "libcurl"
  )
}

download_manifest <- data.frame(
  dataset_id = dataset_id,
  geo_accession = geo_accession,
  source = "GEO",
  source_url = paste0(geo_base_url, geo_files),
  local_file = local_files,
  download_time = as.character(Sys.time())
)

write.csv(
  download_manifest,
  file.path(processed_dir, "download_manifest.csv"),
  row.names = FALSE
)
