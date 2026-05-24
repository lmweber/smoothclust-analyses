################################################################################
# Script to download and extract dataset
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from root directory of smoothclust-analyses repo

dataset_id <- "gbm_greenwald_2024"
geo_accession <- "GSE237183"

repo_dir <- "."
data_dir <- file.path(repo_dir, "data", dataset_id)
raw_dir <- file.path(data_dir, "raw")
extracted_dir <- file.path(data_dir, "extracted")
processed_dir <- file.path(data_dir, "processed")

geo_raw_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/",
  geo_accession,
  "/suppl/",
  geo_accession,
  "_RAW.tar"
)

geo_series_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
  geo_accession
)

# alternative links - not used
# data
zenodo_url <- "https://zenodo.org/records/12624860"
# original analysis code
github_url <- "https://github.com/tiroshlab/Spatial_Glioma"

raw_archive <- file.path(raw_dir, paste0(geo_accession, "_RAW.tar"))


# download data

download.file(
  url = geo_raw_url,
  destfile = raw_archive,
  mode = "wb",
  method = "auto"
)


# extract data

utils::untar(raw_archive, exdir = extracted_dir)


# GEO extracts files in a flat structure. The preprocessing script will load
# directly from these flat files.

extracted_files <- list.files(
  extracted_dir,
  full.names = TRUE,
  recursive = FALSE
)

n_files <- length(extracted_files)

download_manifest <- data.frame(
  dataset_id = dataset_id,
  geo_accession = geo_accession,
  source = "GEO",
  source_url = geo_raw_url,
  local_file = raw_archive,
  extracted_dir = extracted_dir,
  n_extracted_files = n_files,
  download_or_check_time = as.character(Sys.time()),
  stringsAsFactors = FALSE
)

write.csv(
  download_manifest,
  file = file.path(processed_dir, "download_manifest.csv"),
  row.names = FALSE
)
