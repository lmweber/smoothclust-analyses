################################################################################
# Script to evaluate author-derived spot labels
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from root directory of smoothclust-analyses repo

library(smoothclust)
library(SpatialExperiment)

dataset_id <- "GBM"

repo_dir <- "."
data_dir <- file.path(repo_dir, "data", dataset_id)
results_dir <- file.path(repo_dir, "results", dataset_id, "boundary_density")

smoothclust_spe_plot_file <- file.path(results_dir, "spe_GBM_BD_plot.rds")
author_metadata_file <- file.path(data_dir, "raw", "Inputs", "general", "visium_metadata.csv")

# author-derived label sets:
# (i) BD metrics calculated on org1 labels
# (ii) BD metrics calculated on org2 labels
# (iii) BD metrics calculated on Leiden/metaprogram labels
author_label_metadata <- data.frame(
  author_label_col = c("org1", "org2", "mp"),
  author_label_section = c(
    "org1_zonal_labels",
    "org2_zonal_labels",
    "leiden_labels"
  ),
  author_label_display = c(
    "org1 zones (3 classes)",
    "org2 zones (5 classes)",
    "Leiden metaprogram labels"
  ),
  author_label_description = c(
    "structured_disorganized_intermediate_zones",
    "structured_disorganized_malignancy_zones",
    "leiden_metaprogram_label"
  ),
  stringsAsFactors = FALSE
)

author_label_cols <- author_label_metadata$author_label_col
boundary_k <- 6

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# load data

spe <- readRDS(smoothclust_spe_plot_file)
author_metadata <- read.csv(author_metadata_file, stringsAsFactors = FALSE)

stopifnot(all(c("sample", "spot_id", author_label_cols) %in% names(author_metadata)))

sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))


# map sample identifiers between the processed GEO objects and author metadata

author_sample_map <- data.frame(
  sample_id = c(
    "bwh23", "bwh24", "bwh25", "bwh28", "bwh35",
    "mgh258", "mgh259",
    "zh1007inf", "zh1007nec", "zh1019inf", "zh1019t1",
    "zh8811a", "zh8811b", "zh8812", "zh881inf", "zh881t1",
    "zh916bulk", "zh916inf", "zh916t1"
  ),
  author_sample_id = c(
    "bwh23", "bwh24", "bwh25", "bwh28", "bwh35",
    "MGH258", "mgh259",
    "ZH1007inf", "ZH1007nec", "ZH1019inf", "ZH1019T1",
    "ZH8811Abulk", "ZH8811Bbulk", "ZH8812bulk", "ZH881inf", "ZH881T1",
    "ZH916bulk", "ZH916inf", "ZH916T1"
  ),
  stringsAsFactors = FALSE
)

stopifnot(all(sample_ids %in% author_sample_map$sample_id))

sample_match <- match(as.character(colData(spe)$sample_id), author_sample_map$sample_id)
colData(spe)$author_sample_id <- author_sample_map$author_sample_id[sample_match]

author_key <- paste(author_metadata$sample, author_metadata$spot_id, sep = "___")
stopifnot(!anyDuplicated(author_key))

spe_key <- paste(colData(spe)$author_sample_id, colData(spe)$barcode, sep = "___")
author_match <- match(spe_key, author_key)
colData(spe)$author_metadata_match <- !is.na(author_match)

for (author_label_col in author_label_cols) {
  output_col <- paste0("author_", author_label_col)
  colData(spe)[[output_col]] <- author_metadata[[author_label_col]][author_match]
}


# summarize author label coverage

coverage_summary <- data.frame()

for (sample_id in sample_ids) {
  sample_ix <- as.character(colData(spe)$sample_id) == sample_id
  author_sample_id <- unique(as.character(colData(spe)$author_sample_id)[sample_ix])[1]
  
  for (author_label_col in author_label_cols) {
    author_col <- paste0("author_", author_label_col)
    labels <- colData(spe)[[author_col]][sample_ix]
    label_ix <- match(author_label_col, author_label_metadata$author_label_col)
    
    coverage_summary <- rbind(
      coverage_summary,
      data.frame(
        sample_id = sample_id,
        geo_sample_id = unique(as.character(colData(spe)$geo_sample_id)[sample_ix])[1],
        author_sample_id = author_sample_id,
        author_label_col = author_label_col,
        author_label_section = author_label_metadata$author_label_section[label_ix],
        author_label_display = author_label_metadata$author_label_display[label_ix],
        author_label_description = author_label_metadata$author_label_description[label_ix],
        n_spots = sum(sample_ix),
        n_author_metadata_spots = sum(author_metadata$sample == author_sample_id),
        n_labeled_spots = sum(!is.na(labels)),
        prop_labeled_spots = sum(!is.na(labels)) / sum(sample_ix),
        n_author_labels = length(unique(labels[!is.na(labels)])),
        stringsAsFactors = FALSE
      )
    )
  }
}


# calculate boundary density summaries on labeled spots

author_bd_summary <- data.frame()

for (sample_id in sample_ids) {
  message("Processing sample: ", sample_id)
  sample_ix <- as.character(colData(spe)$sample_id) == sample_id
  n_spots <- sum(sample_ix)
  geo_sample_id <- unique(as.character(colData(spe)$geo_sample_id)[sample_ix])[1]
  author_sample_id <- unique(as.character(colData(spe)$author_sample_id)[sample_ix])[1]
  
  for (author_label_col in author_label_cols) {
    label_ix <- match(author_label_col, author_label_metadata$author_label_col)
    author_col <- paste0("author_", author_label_col)
    local_col <- paste0("local_bd_", author_col)
    labeled_ix <- sample_ix & !is.na(colData(spe)[[author_col]])
    
    if (sum(labeled_ix) <= boundary_k) {
      next
    }
    
    spe_sub <- spe[, labeled_ix]
    author_labels <- factor(colData(spe_sub)[[author_col]])
    
    bd <- boundary_density(
      spatialCoords(spe_sub),
      author_labels,
      k = boundary_k,
      adjust = "analytic"
    )
    
    local_values <- rep(NA_real_, ncol(spe))
    if (local_col %in% names(colData(spe))) {
      local_values <- colData(spe)[[local_col]]
    }
    local_values[labeled_ix] <- bd$local_boundary_density
    colData(spe)[[local_col]] <- local_values
    
    author_bd_summary <- rbind(
      author_bd_summary,
      data.frame(
        sample_id = sample_id,
        geo_sample_id = geo_sample_id,
        author_sample_id = author_sample_id,
        method = "author_metadata",
        author_label_col = author_label_col,
        author_label_section = author_label_metadata$author_label_section[label_ix],
        author_label_display = author_label_metadata$author_label_display[label_ix],
        author_label_description = author_label_metadata$author_label_description[label_ix],
        label_col = author_col,
        n_labels = length(unique(author_labels)),
        n_spots = n_spots,
        n_labeled_spots = ncol(spe_sub),
        boundary_k = boundary_k,
        mean_discordant = bd$mean_discordant,
        boundary_density = bd$boundary_density,
        expected_boundary_density = bd$expected_boundary_density,
        relative_boundary_density = bd$relative_boundary_density,
        excess_boundary_density = bd$excess_boundary_density,
        stringsAsFactors = FALSE
      )
    )
  }
}


# save outputs

parameters <- data.frame(
  dataset_id = dataset_id,
  smoothclust_spe_plot_file = smoothclust_spe_plot_file,
  author_metadata_file = author_metadata_file,
  author_label_cols = paste(author_label_cols, collapse = ","),
  author_label_sections = paste(author_label_metadata$author_label_section, collapse = ","),
  author_label_displays = paste(author_label_metadata$author_label_display, collapse = ","),
  comparison_scope = "author_labeled_spots",
  boundary_k = boundary_k,
  stringsAsFactors = FALSE
)

write.csv(
  coverage_summary,
  file = file.path(results_dir, "GBM_BD_author_labels_coverage.csv"),
  row.names = FALSE
)

write.csv(
  author_bd_summary,
  file = file.path(results_dir, "GBM_BD_author_labels_summary.csv"),
  row.names = FALSE
)

write.csv(
  parameters,
  file = file.path(results_dir, "GBM_BD_author_labels_parameters.csv"),
  row.names = FALSE
)


# create and save lightweight SPE object for plotting

spe_plot <- spe[0, ]
reducedDims(spe_plot) <- S4Vectors::SimpleList()

saveRDS(
  spe_plot,
  file = file.path(results_dir, "spe_GBM_author_labels_BD_plot.rds")
)
