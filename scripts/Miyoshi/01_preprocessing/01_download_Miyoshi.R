################################################################################
# Download the Miyoshi et al. human Alzheimer disease Visium dataset
################################################################################

# Run from the smoothclust-analyses repository root directory.

dataset_id <- "Miyoshi"
raw_dir <- file.path("data", dataset_id, "raw")
processed_dir <- file.path("data", dataset_id, "processed")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

geo_base <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233208/suppl/"
)
geo_file <- "GSE233208_Human_visium_ADDS_seurat_processed.rds.gz"

files <- data.frame(
  type = "author_processed_seurat",
  source = "GEO GSE233208",
  url = paste0(geo_base, geo_file),
  file = file.path(raw_dir, geo_file)
)

options(timeout = 7200)
for (i in seq_len(nrow(files))) {
  if (!file.exists(files$file[i])) {
    message("Downloading ", basename(files$file[i]))
    part <- paste0(files$file[i], ".part")
    stopifnot(download.file(
      files$url[i], part, mode = "wb", method = "libcurl"
    ) == 0L)
    stopifnot(file.rename(part, files$file[i]))
  }
}

files$bytes <- file.info(files$file)$size
files$md5 <- unname(tools::md5sum(files$file))
write.csv(
  files, file.path(processed_dir, "download_manifest.csv"), row.names = FALSE
)
