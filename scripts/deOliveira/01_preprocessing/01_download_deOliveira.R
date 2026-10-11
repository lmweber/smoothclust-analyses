################################################################################
# Download one Visium HD section and one Xenium section from de Oliveira et al.
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(BiocFileCache)
library(OSTA.data)

dataset_id <- "deOliveira"
raw_dir <- file.path("data", dataset_id, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(normalizePath(raw_dir), "BiocFileCache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

osta_ids <- c(
  "VisiumHD_HumanColon_Oliveira",
  "Xenium_HumanColon_Oliveira"
)
bfc <- BiocFileCache(cache_dir, ask = FALSE)
repo_dir <- getwd()

for (osta_id in osta_ids) {
  archive <- OSTA.data_load(
    osta_id,
    bfc = bfc,
    pol = FALSE,
    mol = FALSE
  )
  setwd(repo_dir)
  out_dir <- file.path(raw_dir, osta_id)
  dir.create(out_dir, showWarnings = FALSE)
  if (!length(list.files(out_dir, all.files = FALSE))) {
    unzip(archive, exdir = out_dir)
  }
}

downloaded <- unlist(
  lapply(
    file.path(raw_dir, osta_ids),
    list.files,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  ),
  use.names = FALSE
)
manifest <- data.frame(
  file = downloaded,
  bytes = file.info(downloaded)$size,
  md5 = unname(tools::md5sum(downloaded))
)
write.csv(
  manifest,
  file.path(raw_dir, "download_manifest.csv"),
  row.names = FALSE
)
