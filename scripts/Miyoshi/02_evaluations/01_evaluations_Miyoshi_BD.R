################################################################################
# Evaluate boundary density on Miyoshi et al. cortical-region annotations
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(smoothclust)
library(SpatialExperiment)

dataset_id <- "Miyoshi"
boundary_k <- 6L
n_threads <- 4L
input_file <- file.path("data", dataset_id, "processed", "spe_Miyoshi.rds")
results_dir <- file.path("results", dataset_id)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

spe <- readRDS(input_file)
spots <- as.data.frame(colData(spe))
spots$x <- spatialCoords(spe)[, 1]
spots$y <- spatialCoords(spe)[, 2]
schemes <- c("author_cortical_region", "author_network_region")
section_ids <- unique(spots$section_id)
summary_list <- count_list <- local_list <- list()

# Build each graph within one section, using spots with nonmissing labels.
result_i <- 1L
for (scheme in schemes) {
  for (section_id in section_ids) {
    message("Processing ", scheme, ": ", section_id)
    section_ix <- spots$section_id == section_id
    labels <- spots[[scheme]][section_ix]
    keep <- !is.na(labels) & labels != ""
    bd <- boundary_density(
      as.matrix(spots[section_ix, c("x", "y")])[keep, , drop = FALSE],
      labels[keep], k = boundary_k, n_threads = n_threads,
      adjust = "analytic"
    )
    label_counts <- table(labels[keep])
    section_metadata <- spots[which(section_ix)[1], ]

    summary_list[[result_i]] <- data.frame(
      section_id = section_id,
      sample_id = section_metadata$sample_id,
      label_scheme = scheme,
      n_spots = sum(section_ix),
      n_labeled_spots = sum(keep),
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
    count_list[[result_i]] <- data.frame(
      section_id = section_id, label_scheme = scheme,
      label = names(label_counts), n_spots = as.integer(label_counts),
      stringsAsFactors = FALSE
    )
    local_ix <- which(section_ix)[keep]
    local_list[[result_i]] <- data.frame(
      spots[local_ix, c(
        "spot_id", "section_id", "sample_id", "x", "y"
      )],
      label_scheme = scheme,
      author_label = labels[keep],
      local_boundary_density = bd$local_boundary_density,
      stringsAsFactors = FALSE
    )
    result_i <- result_i + 1L
  }
}

bd_summary <- do.call(rbind, summary_list)
label_counts <- do.call(rbind, count_list)
local_spots <- do.call(rbind, local_list)

scheme_summary <- do.call(rbind, lapply(
  split(bd_summary, bd_summary$label_scheme),
  function(x) data.frame(
    label_scheme = x$label_scheme[1],
    n_sections = nrow(x),
    min_relative_boundary_density = min(x$relative_boundary_density),
    q1_relative_boundary_density = unname(quantile(
      x$relative_boundary_density, 0.25
    )),
    median_relative_boundary_density = median(x$relative_boundary_density),
    q3_relative_boundary_density = unname(quantile(
      x$relative_boundary_density, 0.75
    )),
    max_relative_boundary_density = max(x$relative_boundary_density),
    median_boundary_density = median(x$boundary_density),
    stringsAsFactors = FALSE
  )
))
rownames(scheme_summary) <- NULL

# Paired differences are WM-collapsed minus the deposited q=9 annotation.
q9 <- bd_summary[
  bd_summary$label_scheme == "author_cortical_region",
  c("section_id", "boundary_density", "relative_boundary_density")
]
network <- bd_summary[
  bd_summary$label_scheme == "author_network_region",
  c("section_id", "boundary_density", "relative_boundary_density")
]
names(q9)[-1] <- paste0(names(q9)[-1], "_q9")
names(network)[-1] <- paste0(names(network)[-1], "_wm_collapsed")
scheme_comparison <- merge(q9, network, by = "section_id", sort = FALSE)
scheme_comparison$change_relative_boundary_density <-
  scheme_comparison$relative_boundary_density_wm_collapsed -
  scheme_comparison$relative_boundary_density_q9
scheme_comparison$change_boundary_density <-
  scheme_comparison$boundary_density_wm_collapsed -
  scheme_comparison$boundary_density_q9
relative_bd_change <- scheme_comparison$change_relative_boundary_density
comparison_summary <- data.frame(
  n_sections = nrow(scheme_comparison),
  n_decreased = sum(relative_bd_change < 0),
  n_unchanged = sum(relative_bd_change == 0),
  median_change_relative_boundary_density = median(relative_bd_change),
  mean_change_relative_boundary_density = mean(relative_bd_change),
  spearman_rank_correlation = cor(
    scheme_comparison$relative_boundary_density_q9,
    scheme_comparison$relative_boundary_density_wm_collapsed,
    method = "spearman"
  )
)

parameters <- data.frame(
  dataset_id = dataset_id,
  label_scheme = schemes,
  label_interpretation = c(
    "published nine-region BayesSpace q=9 annotation",
    "q=9 annotation with WM1-WM3 collapsed to a single WM label"
  ),
  graph = paste(
    "six-nearest-neighbor graph on the transformed Visium hexagonal lattice;",
    "at tissue edges or gaps, neighbors need not be immediate lattice neighbors"
  ),
  boundary_k = boundary_k,
  adjustment = "analytic",
  deposited_spots = ncol(spe),
  article_reported_spots = 115251L,
  spot_count_note = paste(
    "The deposited object and Supplementary Table 1 contain 115451 spots,",
    "200 more than the 115251 stated in the article."
  ),
  n_threads = n_threads,
  stringsAsFactors = FALSE
)

write.csv(
  bd_summary, file.path(results_dir, "Miyoshi_BD_section_summary.csv"),
  row.names = FALSE
)
write.csv(
  scheme_summary, file.path(results_dir, "Miyoshi_BD_descriptive_summary.csv"),
  row.names = FALSE
)
write.csv(
  scheme_comparison, file.path(results_dir, "Miyoshi_BD_scheme_comparison.csv"),
  row.names = FALSE
)
write.csv(
  comparison_summary,
  file.path(results_dir, "Miyoshi_BD_scheme_comparison_summary.csv"),
  row.names = FALSE
)
write.csv(
  label_counts, file.path(results_dir, "Miyoshi_author_label_counts.csv"),
  row.names = FALSE
)
saveRDS(local_spots, file.path(results_dir, "Miyoshi_BD_spot_data.rds"))
write.csv(
  parameters, file.path(results_dir, "Miyoshi_BD_parameters.csv"),
  row.names = FALSE
)
