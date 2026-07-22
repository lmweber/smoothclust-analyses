################################################################################
# Download the Kuppe et al. human myocardial infarction Visium dataset
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(jsonlite)

dataset_id <- "Kuppe"
raw_dir <- file.path("data", dataset_id, "raw")
processed_dir <- file.path("data", dataset_id, "processed")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

supplement_url <- paste0(
  "https://static-content.springer.com/esm/",
  "art%3A10.1038%2Fs41586-022-05060-x/MediaObjects/"
)
files <- data.frame(
  type = c("clinical_metadata", "visium_metadata"),
  source = "Nature supplementary tables",
  url = paste0(
    supplement_url,
    c(
      "41586_2022_5060_MOESM6_ESM.xlsx",
      "41586_2022_5060_MOESM9_ESM.xlsx"
    )
  ),
  file = file.path(
    raw_dir,
    c(
      "Kuppe2022_Supplementary_Table_1.xlsx",
      "Kuppe2022_Supplementary_Table_4.xlsx"
    )
  )
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

# The contributor Seurat objects retain the authors' exact cell-type-niche and
# molecular-niche spot labels, together with counts and spatial coordinates.
hca_project <- "e9f36305-d857-44a3-93f0-df4e6007dc97"
hca_filters <- URLencode(paste0(
  '{"projectId":{"is":["', hca_project,
  '"]},"fileFormat":{"is":["rds"]}}'
), reserved = TRUE)
hca_index_url <- paste0(
  "https://service.azul.data.humancellatlas.org/index/files",
  "?catalog=dcp60&size=100&filters=", hca_filters
)
hca <- fromJSON(hca_index_url, simplifyVector = FALSE)
rds_entries <- lapply(hca$hits, function(x) x$files[[1]])
rds_dir <- file.path(raw_dir, "hca_rds")
dir.create(rds_dir, showWarnings = FALSE)
rds_files <- data.frame(
  type = "author_annotated_visium_rds",
  source = paste("Human Cell Atlas project", hca_project),
  url = vapply(rds_entries, `[[`, character(1), "azul_url"),
  file = file.path(
    rds_dir, vapply(rds_entries, `[[`, character(1), "name")
  )
)
files <- rbind(files, rds_files)

for (i in seq.int(3L, nrow(files))) {
  if (!file.exists(files$file[i])) {
    message("Downloading ", basename(files$file[i]), " (", i - 2L, "/",
            nrow(files) - 2L, ")")
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
