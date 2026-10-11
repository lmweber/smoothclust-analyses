################################################################################
# Evaluate boundary density on Maynard et al. manual layer annotations
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(smoothclust)

dataset_id <- "Maynard"
boundary_k <- 6L
n_threads <- 4L
input_file <- file.path(
  "data", dataset_id, "processed", "Maynard_spot_metadata.rds"
)
results_dir <- file.path("results", dataset_id)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

spots <- readRDS(input_file)
section_ids <- unique(spots$section_id)
summary_list <- count_list <- local_list <- vector("list", length(section_ids))

for (i in seq_along(section_ids)) {
  section_id <- section_ids[i]
  message("Processing ", section_id)
  section_ix <- spots$section_id == section_id
  labels <- spots$author_layer[section_ix]
  # Exclude missing or empty labels; retain annotated spots marked discard=TRUE.
  keep <- !is.na(labels) & labels != ""
  bd <- boundary_density(
    as.matrix(spots[section_ix, c("x", "y")])[keep, , drop = FALSE],
    labels[keep], k = boundary_k, n_threads = n_threads,
    adjust = "analytic"
  )
  label_counts <- table(labels[keep])
  section_metadata <- spots[which(section_ix)[1], ]

  summary_list[[i]] <- data.frame(
    section_id = section_id,
    sample_id = section_metadata$sample_id,
    patient_id = section_metadata$patient_id,
    position_um = section_metadata$position_um,
    adjacent_replicate = section_metadata$adjacent_replicate,
    pair_id = section_metadata$pair_id,
    n_spots = sum(section_ix),
    n_labeled_spots = sum(keep),
    n_author_discard = sum(spots$author_discard[section_ix], na.rm = TRUE),
    n_labeled_author_discard = sum(
      keep & spots$author_discard[section_ix], na.rm = TRUE
    ),
    labeled_fraction = mean(keep),
    n_labels = length(label_counts),
    min_label_size = min(label_counts),
    boundary_k = boundary_k,
    mean_discordant = bd$mean_discordant,
    boundary_density = bd$boundary_density,
    expected_boundary_density = bd$expected_boundary_density,
    relative_boundary_density = bd$relative_boundary_density,
    stringsAsFactors = FALSE
  )
  count_list[[i]] <- data.frame(
    section_id = section_id,
    label = names(label_counts),
    n_spots = as.integer(label_counts),
    stringsAsFactors = FALSE
  )
  local_ix <- which(section_ix)[keep]
  local_list[[i]] <- data.frame(
    spots[local_ix, c(
      "spot_id", "section_id", "sample_id", "patient_id", "position_um",
      "adjacent_replicate", "pair_id", "x", "y", "author_layer"
    )],
    local_boundary_density = bd$local_boundary_density,
    stringsAsFactors = FALSE
  )
}

bd_summary <- do.call(rbind, summary_list)
label_counts <- do.call(rbind, count_list)
local_spots <- do.call(rbind, local_list)

# Pair replicates 1 and 2 within each donor and tissue-block position.
pair_wide <- reshape(
  bd_summary[c(
    "pair_id", "patient_id", "position_um", "adjacent_replicate",
    "section_id", "boundary_density", "relative_boundary_density"
  )],
  idvar = c("pair_id", "patient_id", "position_um"),
  timevar = "adjacent_replicate", direction = "wide"
)
pair_wide$absolute_difference_boundary_density <- abs(
  pair_wide$boundary_density.2 - pair_wide$boundary_density.1
)
pair_wide$absolute_difference_relative_boundary_density <- abs(
  pair_wide$relative_boundary_density.2 -
    pair_wide$relative_boundary_density.1
)

descriptive_summary <- data.frame(
  n_sections = nrow(bd_summary),
  n_patients = length(unique(bd_summary$patient_id)),
  n_adjacent_pairs = nrow(pair_wide),
  min_relative_boundary_density = min(bd_summary$relative_boundary_density),
  median_relative_boundary_density = median(
    bd_summary$relative_boundary_density
  ),
  max_relative_boundary_density = max(bd_summary$relative_boundary_density),
  median_boundary_density = median(bd_summary$boundary_density),
  median_adjacent_pair_absolute_difference = median(
    pair_wide$absolute_difference_relative_boundary_density
  ),
  max_adjacent_pair_absolute_difference = max(
    pair_wide$absolute_difference_relative_boundary_density
  )
)
parameters <- data.frame(
  dataset_id = dataset_id,
  label = "published manual six-layer plus white-matter annotation",
  graph = paste(
    "six-nearest-neighbor graph on the transformed Visium hexagonal lattice;",
    "at tissue edges or gaps, neighbors need not be immediate lattice neighbors"
  ),
  boundary_k = boundary_k,
  adjustment = "analytic",
  missing_labels = "352 author-unassigned spots excluded",
  discard_handling = paste(
    "534 spots with a nonmissing manual layer label and deposited",
    "discard=TRUE were retained; the analysis summarizes the geometry",
    "of the provided manual labels"
  ),
  n_threads = n_threads,
  stringsAsFactors = FALSE
)

write.csv(
  bd_summary, file.path(results_dir, "Maynard_BD_section_summary.csv"),
  row.names = FALSE
)
write.csv(
  pair_wide, file.path(results_dir, "Maynard_BD_adjacent_pair_summary.csv"),
  row.names = FALSE
)
write.csv(
  descriptive_summary,
  file.path(results_dir, "Maynard_BD_descriptive_summary.csv"),
  row.names = FALSE
)
write.csv(
  label_counts, file.path(results_dir, "Maynard_author_label_counts.csv"),
  row.names = FALSE
)
saveRDS(local_spots, file.path(results_dir, "Maynard_BD_spot_data.rds"))
write.csv(
  parameters, file.path(results_dir, "Maynard_BD_parameters.csv"),
  row.names = FALSE
)
