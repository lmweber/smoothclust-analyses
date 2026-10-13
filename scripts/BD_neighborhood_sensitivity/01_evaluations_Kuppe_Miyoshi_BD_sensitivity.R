################################################################################
# BD neighborhood sensitivity for fixed Kuppe and Miyoshi published labels
################################################################################

# Run from the smoothclust-analyses repository root directory.
# Read processed inputs; do not rerun preprocessing, smoothing, or clustering.

library(smoothclust)
library(SpatialExperiment)

boundary_k_values <- c(6L, 18L, 36L)
n_threads <- 4L
# Absolute tolerance accommodates the saved CSVs' decimal serialization.
baseline_tolerance <- 1e-12
results_dir <- file.path("results", "BD_neighborhood_sensitivity")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- data.frame(
  dataset_id = c("Kuppe", "Miyoshi"),
  label_scheme = c("author_celltype_niche", "author_cortical_region"),
  label_definition = c(
    "Nine published computational cell-type-niche labels",
    "Original nine-region BayesSpace q=9 annotation"
  ),
  baseline_prefix = c("Kuppe_celltype_niche_BD", "Miyoshi_BD"),
  counts_file = c(
    "Kuppe_author_celltype_niche_counts.csv", "Miyoshi_author_label_counts.csv"
  ),
  expected_sections = c(28L, 39L),
  stringsAsFactors = FALSE
)
metadata_columns <- c(
  "sample_id", "patient_id", "sample_group", "author_group_code",
  "region_code", "region_label", "region_novel", "days_after_infarction",
  "disease", "age", "sex", "author_image_id"
)
metrics <- c(
  "boundary_density", "expected_boundary_density", "relative_boundary_density"
)
summaries <- list()

same_values <- function(x, y) {
  length(x) == length(y) && all(
    (is.na(x) & is.na(y)) |
      (is.finite(x) & is.finite(y) & abs(x - y) <= baseline_tolerance)
  )
}

for (dataset_i in seq_len(nrow(datasets))) {
  config <- datasets[dataset_i, ]
  dataset_id <- config$dataset_id
  scheme <- config$label_scheme
  baseline_dir <- file.path("results", dataset_id)
  input_files <- c(
    processed = file.path("data", dataset_id, "processed",
                          paste0("spe_", dataset_id, ".rds")),
    summary = file.path(baseline_dir,
                        paste0(config$baseline_prefix, "_section_summary.csv")),
    spots = file.path(baseline_dir,
                      paste0(config$baseline_prefix, "_spot_data.rds")),
    counts = file.path(baseline_dir, config$counts_file)
  )
  if (any(!file.exists(input_files))) {
    stop("Missing inputs: ", paste(input_files[!file.exists(input_files)],
                                   collapse = "; "), call. = FALSE)
  }
  spe <- readRDS(input_files[["processed"]])
  spots <- as.data.frame(colData(spe))
  spots$x <- spatialCoords(spe)[, 1]
  spots$y <- spatialCoords(spe)[, 2]
  baseline <- read.csv(input_files[["summary"]], stringsAsFactors = FALSE)
  baseline <- baseline[baseline$label_scheme == scheme, ]
  saved_spots <- readRDS(input_files[["spots"]])
  saved_counts <- read.csv(input_files[["counts"]], stringsAsFactors = FALSE)
  if (dataset_id == "Kuppe") {
    saved_spots$author_label <- saved_spots[[scheme]]
    names(saved_counts)[names(saved_counts) == scheme] <- "label"
  } else {
    saved_spots <- saved_spots[saved_spots$label_scheme == scheme, ]
    saved_counts <- saved_counts[saved_counts$label_scheme == scheme, ]
  }
  retained <- !is.na(saved_spots$author_label) & saved_spots$author_label != ""
  saved_spots <- saved_spots[retained, ]
  section_ids <- unique(as.character(spots$section_id))
  stopifnot(
    "Section identities must match the saved baseline" =
      length(section_ids) == config$expected_sections &&
      identical(section_ids, as.character(baseline$section_id)),
    "Saved baseline must use k=6" = all(baseline$boundary_k == 6L)
  )
  if (dataset_id == "Kuppe") {
    stopifnot("Kuppe categories must be present" =
                !anyNA(spots$sample_group) && all(spots$sample_group != ""))
  }

  for (section_id in section_ids) {
    message("Processing ", dataset_id, ": ", section_id)
    section <- spots[spots$section_id == section_id, ]
    labels <- section[[scheme]]
    keep <- !is.na(labels) & labels != ""
    current <- section[keep, ]
    labels <- labels[keep]
    coords <- as.matrix(current[, c("x", "y")])
    old <- saved_spots[saved_spots$section_id == section_id, ]
    old_summary <- baseline[baseline$section_id == section_id, ]
    old_counts <- saved_counts[saved_counts$section_id == section_id, ]
    label_counts <- table(labels)
    stopifnot(
      "Retained spot identities, order, labels and coordinates must match" =
        identical(as.character(current$spot_id), as.character(old$spot_id)) &&
        identical(as.character(labels), as.character(old$author_label)) &&
        identical(unname(coords), unname(as.matrix(old[, c("x", "y")]))),
      "Label counts must match" =
        identical(names(label_counts), as.character(old_counts$label)) &&
        identical(as.integer(label_counts), as.integer(old_counts$n_spots)),
      "Each section must retain more than 36 spots" =
        nrow(current) > max(boundary_k_values)
    )

    row <- data.frame(
      dataset_id = dataset_id, section_id = section_id,
      label_scheme = scheme, label_definition = config$label_definition,
      n_spots = nrow(section), n_labeled_spots = nrow(current),
      n_labels = length(label_counts), stringsAsFactors = FALSE
    )
    for (column in metadata_columns) {
      row[[column]] <- if (column %in% names(section)) section[[column]][1] else NA
    }
    for (boundary_k in boundary_k_values) {
      # Retain input order and the baseline's default search conventions.
      # Baseline checks expose differences in ties across package environments.
      # findKNN excludes the focal location; do not truncate a maximum-k graph.
      bd <- boundary_density(coords, labels, k = boundary_k,
                             n_threads = n_threads, adjust = TRUE)
      if (boundary_k == 6L) {
        for (metric in metrics) {
          if (!same_values(bd[[metric]], old_summary[[metric]])) {
            stop(dataset_id, " / ", section_id, ": k=6 ", metric,
                 " differs from baseline (recomputed: ", bd[[metric]],
                 "; saved: ", old_summary[[metric]], "; tolerance: ",
                 baseline_tolerance, ")", call. = FALSE)
          }
        }
        stopifnot("k=6 local BD must match the saved baseline" =
                    same_values(bd$local_boundary_density,
                                old$local_boundary_density))
        baseline_expected <- bd$expected_boundary_density
      } else {
        stopifnot("Expected BD must remain unchanged across k" =
                    same_values(bd$expected_boundary_density, baseline_expected))
      }
      result <- row
      result$boundary_k <- boundary_k
      result$mean_discordant <- bd$mean_discordant
      for (metric in metrics) result[[metric]] <- bd[[metric]]
      summaries[[length(summaries) + 1L]] <- result
    }
  }
  rm(spe, spots, saved_spots)
}
section_summary <- do.call(rbind, summaries)
rownames(section_summary) <- NULL

# Undefined values remain in the tables and invalidate affected summaries.
undefined <- section_summary[
  !is.finite(section_summary$boundary_density) |
    !is.finite(section_summary$expected_boundary_density) |
    !is.finite(section_summary$relative_boundary_density),
  c("dataset_id", "section_id", "boundary_k", metrics)
]
undefined$reason <- ifelse(
  is.finite(undefined$expected_boundary_density) &
    undefined$expected_boundary_density == 0,
  "Relative BD undefined: zero fixed-count expectation",
  "Nonfinite BD output; inspect inputs and boundary_density"
)
if (nrow(undefined)) {
  warning("Undefined BD values; see BD_k_sensitivity_undefined.csv:\n",
          paste(undefined$dataset_id, undefined$section_id,
                paste0("k=", undefined$boundary_k), collapse = "\n"),
          call. = FALSE)
}

changes <- comparisons <- list()
for (dataset_id in datasets$dataset_id) {
  base <- section_summary[section_summary$dataset_id == dataset_id &
                            section_summary$boundary_k == 6L, ]
  for (boundary_k in boundary_k_values[-1]) {
    larger <- section_summary[section_summary$dataset_id == dataset_id &
                                section_summary$boundary_k == boundary_k, ]
    change <- larger[, c("dataset_id", "section_id", metadata_columns,
                         "label_scheme", "boundary_k")]
    for (metric in c("boundary_density", "relative_boundary_density")) {
      delta <- larger[[metric]] - base[[metric]]
      defined <- is.finite(base[[metric]]) & is.finite(larger[[metric]])
      change[[paste0(metric, "_k6")]] <- base[[metric]]
      change[[paste0(metric, "_larger_k")]] <- larger[[metric]]
      change[[paste0("change_", metric)]] <- delta
      complete <- all(defined)
      variable <- complete && length(unique(base[[metric]])) > 1L &&
        length(unique(larger[[metric]])) > 1L
      comparisons[[length(comparisons) + 1L]] <- data.frame(
        dataset_id = dataset_id, boundary_k = boundary_k, metric = metric,
        n_sections = nrow(base), n_defined = sum(defined),
        spearman_rank_correlation = if (variable)
          cor(base[[metric]], larger[[metric]], method = "spearman") else NA_real_,
        median_signed_change = if (complete) median(delta) else NA_real_,
        median_absolute_change = if (complete) median(abs(delta)) else NA_real_,
        maximum_absolute_change = if (complete) max(abs(delta)) else NA_real_,
        status = if (!complete) "Undefined values; no sections excluded" else
          if (!variable) "Spearman undefined: constant values" else "Defined",
        stringsAsFactors = FALSE
      )
    }
    changes[[length(changes) + 1L]] <- change
  }
}

# All sections contribute equally; retain dominant-niche sections at every k.
kuppe <- section_summary[section_summary$dataset_id == "Kuppe", ]
category_summary <- do.call(rbind, lapply(
  split(kuppe, list(kuppe$sample_group, kuppe$boundary_k), drop = TRUE),
  function(x) do.call(rbind, lapply(metrics, function(metric) {
    values <- x[[metric]]
    defined <- is.finite(values)
    quartiles <- if (all(defined)) {
      unname(quantile(values, c(0.25, 0.5, 0.75)))
    } else rep(NA_real_, 3)
    data.frame(
      dataset_id = "Kuppe", sample_group = x$sample_group[1],
      boundary_k = x$boundary_k[1], metric = metric,
      n_sections = nrow(x), n_patients = length(unique(x$patient_id)),
      n_defined = sum(defined), q1 = quartiles[1], median = quartiles[2],
      q3 = quartiles[3], iqr = quartiles[3] - quartiles[1],
      stringsAsFactors = FALSE
    )
  }))
))
parameters <- data.frame(
  boundary_k_values = paste(boundary_k_values, collapse = ","),
  n_threads = n_threads, baseline_absolute_tolerance = baseline_tolerance,
  adjustment = "analytic", coordinates = "Unchanged processed spatialCoords",
  graph = "Directed Euclidean findKNN, excluding focal locations",
  ties = "Unchanged input order and independent default findKNN calls at each k",
  scope = "Fixed published labels; descriptive equal-section summaries",
  undefined_values = "Retained; affected summaries are NA without exclusions",
  R_version = R.version.string,
  smoothclust_version = as.character(packageVersion("smoothclust")),
  BiocNeighbors_version = as.character(packageVersion("BiocNeighbors"))
)

# Separate sensitivity files; primary k=6 outputs are never overwritten.
outputs <- list(
  section_summary = section_summary,
  section_changes = do.call(rbind, changes),
  comparison_summary = do.call(rbind, comparisons),
  Kuppe_category_summary = category_summary,
  parameters = parameters
)
for (name in names(outputs)) {
  stem <- file.path(results_dir, paste0("BD_k_sensitivity_", name))
  write.csv(outputs[[name]], paste0(stem, ".csv"), row.names = FALSE)
  saveRDS(outputs[[name]], paste0(stem, ".rds"))
}
write.csv(undefined, file.path(results_dir, "BD_k_sensitivity_undefined.csv"),
          row.names = FALSE)
