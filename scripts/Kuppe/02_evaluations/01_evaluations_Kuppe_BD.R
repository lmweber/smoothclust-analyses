################################################################################
# Evaluate boundary density on Kuppe et al. published cell-type-niche labels
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(smoothclust)
library(SpatialExperiment)

dataset_id <- "Kuppe"
boundary_k <- 6L
n_threads <- 4L
input_file <- file.path("data", dataset_id, "processed", "spe_Kuppe.rds")
results_dir <- file.path("results", dataset_id)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

spe <- readRDS(input_file)
spots <- as.data.frame(colData(spe))
spots$x <- spatialCoords(spe)[, 1]
spots$y <- spatialCoords(spe)[, 2]
spots$local_boundary_density <- NA_real_

section_ids <- unique(spots$section_id)
summary_list <- label_count_list <- vector("list", length(section_ids))

for (i in seq_along(section_ids)) {
  section_id <- section_ids[i]
  message("Processing section ", i, "/", length(section_ids), ": ",
          section_id)
  section_ix <- spots$section_id == section_id
  labels <- spots$author_celltype_niche[section_ix]
  keep <- !is.na(labels) & labels != ""

  bd <- boundary_density(
    as.matrix(spots[section_ix, c("x", "y")])[keep, , drop = FALSE],
    labels[keep],
    k = boundary_k,
    n_threads = n_threads,
    adjust = "analytic"
  )

  local_ix <- which(section_ix)[keep]
  spots$local_boundary_density[local_ix] <- bd$local_boundary_density
  label_counts <- table(labels[keep])
  dominant_niche <- names(label_counts)[which.max(label_counts)]
  dominant_niche_fraction <- max(label_counts) / sum(label_counts)
  section_metadata <- spots[which(section_ix)[1], ]

  summary_list[[i]] <- data.frame(
    section_id = section_id,
    sample_id = section_metadata$sample_id,
    patient_id = section_metadata$patient_id,
    sample_group = section_metadata$sample_group,
    author_group_code = section_metadata$author_group_code,
    region_code = section_metadata$region_code,
    region_label = section_metadata$region_label,
    region_novel = section_metadata$region_novel,
    days_after_infarction = section_metadata$days_after_infarction,
    disease = section_metadata$disease,
    age = section_metadata$age,
    sex = section_metadata$sex,
    label_scheme = "author_celltype_niche",
    n_spots = sum(section_ix),
    n_labeled_spots = sum(keep),
    labeled_fraction = mean(keep),
    n_labels = length(label_counts),
    min_label_size = min(label_counts),
    dominant_niche = dominant_niche,
    dominant_niche_fraction = dominant_niche_fraction,
    dominant_niche_gt80 = dominant_niche_fraction > 0.8,
    boundary_k = boundary_k,
    mean_discordant = bd$mean_discordant,
    boundary_density = bd$boundary_density,
    expected_boundary_density = bd$expected_boundary_density,
    relative_boundary_density = bd$relative_boundary_density,
    stringsAsFactors = FALSE
  )

  label_count_list[[i]] <- data.frame(
    section_id = section_id,
    sample_id = section_metadata$sample_id,
    patient_id = section_metadata$patient_id,
    sample_group = section_metadata$sample_group,
    region_label = section_metadata$region_label,
    author_celltype_niche = names(label_counts),
    n_spots = as.integer(label_counts),
    stringsAsFactors = FALSE
  )
}

bd_summary <- do.call(rbind, summary_list)
label_counts <- do.call(rbind, label_count_list)

summarize_groups <- function(column, grouping) {
  do.call(rbind, lapply(split(bd_summary, bd_summary[[column]]), function(x) {
    data.frame(
      grouping = grouping,
      group = x[[column]][1],
      n_sections = nrow(x),
      n_samples = length(unique(x$sample_id)),
      n_patients = length(unique(x$patient_id)),
      median_boundary_density = median(x$boundary_density),
      q1_boundary_density = unname(quantile(x$boundary_density, 0.25)),
      q3_boundary_density = unname(quantile(x$boundary_density, 0.75)),
      median_expected_boundary_density = median(
        x$expected_boundary_density
      ),
      median_relative_boundary_density = median(x$relative_boundary_density),
      q1_relative_boundary_density = unname(quantile(
        x$relative_boundary_density, 0.25
      )),
      q3_relative_boundary_density = unname(quantile(
        x$relative_boundary_density, 0.75
      )),
      n_sections_dominant_niche_gt80 = sum(x$dominant_niche_gt80),
      dominant_niche_gt80_section_ids = paste(
        x$section_id[x$dominant_niche_gt80], collapse = "; "
      ),
      median_relative_boundary_density_excluding_dominant_niche_gt80 =
        median(x$relative_boundary_density[!x$dominant_niche_gt80]),
      stringsAsFactors = FALSE
    )
  }))
}
group_summary <- rbind(
  summarize_groups(
    "sample_group", "author histomorphological sample category"
  ),
  summarize_groups("region_label", "histomorphological region")
)
rownames(group_summary) <- NULL

parameters <- data.frame(
  dataset_id = dataset_id,
  input_file = input_file,
  label_scheme = "HCA RDS meta.data/celltype_niche",
  label_interpretation = paste(
    "published computational labels: nine cell-type niches derived from",
    "cell2location proportions and joint clustering"
  ),
  sample_group_interpretation = paste(
    "authors' histomorphological sample categories:",
    "myogenic-, ischaemic-, and fibrotic-enriched"
  ),
  coordinates = paste(
    "Visium array positions transformed to an equilateral hexagonal lattice"
  ),
  graph_limitation = paste(
    "a fixed six-nearest-neighbor graph can bridge tissue edges or gaps;",
    "neighbor-count sensitivity was not evaluated"
  ),
  boundary_k = boundary_k,
  adjustment = "analytic",
  dominant_niche_diagnostic = paste(
    "sections with more than 80% of labeled spots assigned to one niche",
    "are flagged; descriptive summaries are also reported without them"
  ),
  missing_labels = "excluded",
  n_threads = n_threads,
  stringsAsFactors = FALSE
)

write.csv(
  bd_summary,
  file.path(results_dir, "Kuppe_celltype_niche_BD_section_summary.csv"),
  row.names = FALSE
)
write.csv(
  group_summary,
  file.path(results_dir, "Kuppe_celltype_niche_BD_group_summary.csv"),
  row.names = FALSE
)
write.csv(
  label_counts,
  file.path(results_dir, "Kuppe_author_celltype_niche_counts.csv"),
  row.names = FALSE
)
saveRDS(
  spots, file.path(results_dir, "Kuppe_celltype_niche_BD_spot_data.rds")
)
write.csv(
  parameters,
  file.path(results_dir, "Kuppe_celltype_niche_BD_parameters.csv"),
  row.names = FALSE
)
