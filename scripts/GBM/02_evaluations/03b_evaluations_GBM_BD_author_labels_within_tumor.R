################################################################################
# Script to evaluate within-tumor author-derived label results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(BiocNeighbors)
library(SpatialExperiment)

dataset_id <- "GBM"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)

spe_plot_file <- file.path(results_dir, "spe_GBM_author_labels_BD_plot.rds")
coverage_file <- file.path(results_dir, "GBM_BD_author_labels_coverage.csv")
author_bd_file <- file.path(results_dir, "GBM_BD_author_labels_summary.csv")
sample_biology_file <- file.path(
  results_dir,
  "GBM_BD_author_labels_sample_biology.csv"
)

section_summary_file <- file.path(
  results_dir,
  "GBM_BD_author_labels_within_tumor_section_summary.csv"
)
composition_file <- file.path(
  results_dir,
  "GBM_BD_author_labels_within_tumor_composition.csv"
)
tumor_summary_file <- file.path(
  results_dir,
  "GBM_BD_author_labels_within_tumor_tumor_summary.csv"
)
qc_file <- file.path(
  results_dir,
  "GBM_BD_author_labels_within_tumor_qc.csv"
)
parameters_file <- file.path(
  results_dir,
  "GBM_BD_author_labels_within_tumor_parameters.csv"
)

boundary_k <- 6L
bd_label_cols <- c("mp", "layer")
coverage_label_cols <- c("org1", "mp", "layer")
validation_tolerance <- 1e-10

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


################################################################################
# Analysis sections
################################################################################

section_metadata <- data.frame(
  sample_id = c(
    "zh1007inf", "zh1007nec",
    "zh1019inf", "zh1019t1",
    "zh8811a", "zh8811b", "zh8812", "zh881inf", "zh881t1",
    "zh916bulk", "zh916inf", "zh916t1"
  ),
  author_sample_id = c(
    "ZH1007inf", "ZH1007nec",
    "ZH1019inf", "ZH1019T1",
    "ZH8811Abulk", "ZH8811Bbulk", "ZH8812bulk", "ZH881inf", "ZH881T1",
    "ZH916bulk", "ZH916inf", "ZH916T1"
  ),
  tumor_group = c(
    "ZH1007", "ZH1007",
    "ZH1019", "ZH1019",
    "ZH881", "ZH881", "ZH881", "ZH881", "ZH881",
    "ZH916", "ZH916", "ZH916"
  ),
  region_group = c(
    "Infiltrating", "Necrotic",
    "Infiltrating", "T1 contrast-enhancing",
    "Bulk tumor", "Bulk tumor", "Bulk tumor", "Infiltrating",
    "T1 contrast-enhancing",
    "Bulk tumor", "Infiltrating", "T1 contrast-enhancing"
  ),
  section_display = c(
    "Infiltrating", "Necrotic",
    "Infiltrating", "T1 contrast-enhancing",
    "1A bulk", "1B bulk", "2 bulk", "Infiltrating",
    "T1 contrast-enhancing",
    "Bulk tumor", "Infiltrating", "T1 contrast-enhancing"
  ),
  tumor_section_order = c(
    1L, 2L,
    1L, 2L,
    1L, 2L, 3L, 4L, 5L,
    1L, 2L, 3L
  ),
  tumor_n_sections = c(
    2L, 2L,
    2L, 2L,
    5L, 5L, 5L, 5L, 5L,
    3L, 3L, 3L
  ),
  stringsAsFactors = FALSE
)

analysis_sample_ids <- section_metadata$sample_id
analysis_tumor_groups <- unique(section_metadata$tumor_group)

if (nrow(section_metadata) != 12L ||
    anyDuplicated(section_metadata$sample_id) ||
    length(analysis_tumor_groups) != 4L) {
  stop("Expected 12 unique sections from four multi-region tumors.")
}


################################################################################
# Helper functions
################################################################################

require_columns <- function(x, required, object_name) {
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0L) {
    stop(
      object_name,
      " is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

analytic_expected_bd <- function(labels) {
  labels <- labels[!is.na(labels)]
  n <- length(labels)
  if (n < 2L) {
    return(NA_real_)
  }
  label_counts <- as.numeric(table(labels))
  1 - sum(label_counts * (label_counts - 1)) / (n * (n - 1))
}

make_directed_knn_edges <- function(spatial_coords, k) {
  if (nrow(spatial_coords) <= k) {
    stop("The number of graph vertices must be greater than boundary_k.")
  }
  nn <- BiocNeighbors::findKNN(
    spatial_coords,
    k = k,
    get.index = TRUE,
    get.distance = TRUE
  )
  data.frame(
    source = rep(seq_len(nrow(nn$index)), times = ncol(nn$index)),
    target = as.vector(nn$index),
    edge_distance = as.vector(nn$distance)
  )
}

summarize_labeled_edges <- function(edges, labels) {
  is_labeled <- !is.na(labels)
  retained <- edges[
    is_labeled[edges$source] & is_labeled[edges$target],
    ,
    drop = FALSE
  ]
  if (nrow(retained) == 0L) {
    stop("No graph edges have two labeled endpoints.")
  }

  retained$discordant <-
    labels[retained$source] != labels[retained$target]
  n_neighbors <- tabulate(retained$source, nbins = length(labels))
  n_discordant <- tabulate(
    retained$source[retained$discordant],
    nbins = length(labels)
  )
  labeled_sources <- which(is_labeled)
  eligible_sources <- labeled_sources[n_neighbors[labeled_sources] > 0L]
  if (length(eligible_sources) == 0L) {
    stop("No labeled source spots retain a labeled neighbor.")
  }

  local_bd <-
    n_discordant[eligible_sources] / n_neighbors[eligible_sources]
  boundary_density_value <- mean(local_bd)
  expected_boundary_density_value <- analytic_expected_bd(labels)
  relative_boundary_density_value <- if (
    is.finite(expected_boundary_density_value) &&
      expected_boundary_density_value > 0
  ) {
    boundary_density_value / expected_boundary_density_value
  } else {
    NA_real_
  }
  retained_neighbors <- n_neighbors[labeled_sources]

  list(
    n_directed_edges_before_filter = nrow(edges),
    n_directed_edges_retained = nrow(retained),
    prop_directed_edges_retained = nrow(retained) / nrow(edges),
    n_labeled_sources_with_neighbors = length(eligible_sources),
    n_labeled_sources_without_neighbors =
      length(labeled_sources) - length(eligible_sources),
    min_retained_neighbors_per_labeled_source = min(retained_neighbors),
    mean_retained_neighbors_per_labeled_source = mean(retained_neighbors),
    median_retained_neighbors_per_labeled_source = median(retained_neighbors),
    boundary_density = boundary_density_value,
    expected_boundary_density = expected_boundary_density_value,
    relative_boundary_density = relative_boundary_density_value,
    mean_edge_distance = mean(retained$edge_distance),
    median_edge_distance = median(retained$edge_distance),
    max_edge_distance = max(retained$edge_distance)
  )
}

make_qc_row <- function(
    section_row,
    author_label_col,
    graph_definition,
    n_graph_vertices,
    graph_summary) {
  data.frame(
    sample_id = section_row$sample_id,
    geo_sample_id = section_row$geo_sample_id,
    author_sample_id = section_row$author_sample_id,
    tumor_group = section_row$tumor_group,
    region_group = section_row$region_group,
    section_display = section_row$section_display,
    tumor_section_order = section_row$tumor_section_order,
    author_label_col = author_label_col,
    graph_definition = graph_definition,
    boundary_k = boundary_k,
    n_section_spots = section_row$n_spots,
    n_graph_vertices = n_graph_vertices,
    n_labeled_spots = section_row[[
      paste0("n_", author_label_col, "_labeled_spots")
    ]],
    prop_labeled_spots = section_row[[
      paste0("prop_", author_label_col, "_labeled_spots")
    ]],
    as.data.frame(graph_summary, stringsAsFactors = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


################################################################################
# Load canonical author-label outputs
################################################################################

input_files <- c(spe_plot_file, coverage_file, author_bd_file, sample_biology_file)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "))
}

spe <- readRDS(spe_plot_file)
coverage <- read.csv(coverage_file, stringsAsFactors = FALSE)
author_bd <- read.csv(author_bd_file, stringsAsFactors = FALSE)
sample_biology <- read.csv(sample_biology_file, stringsAsFactors = FALSE)

require_columns(
  as.data.frame(colData(spe)),
  c("sample_id", "author_mp", "author_layer"),
  "spe colData"
)
require_columns(
  coverage,
  c(
    "sample_id", "geo_sample_id", "author_sample_id", "author_label_col",
    "n_spots", "n_labeled_spots", "prop_labeled_spots", "n_author_labels"
  ),
  "coverage"
)
require_columns(
  author_bd,
  c(
    "sample_id", "author_label_col", "n_labels", "boundary_k",
    "mean_discordant", "boundary_density", "expected_boundary_density",
    "relative_boundary_density"
  ),
  "author_bd"
)
require_columns(
  sample_biology,
  c("sample_id", "fraction_structured", "fraction_mes_hyp"),
  "sample_biology"
)


################################################################################
# Section-level analysis table
################################################################################

section_summary <- section_metadata

for (author_label_col in coverage_label_cols) {
  coverage_sub <- coverage[
    coverage$author_label_col == author_label_col,
    ,
    drop = FALSE
  ]
  coverage_sub <- coverage_sub[
    match(section_summary$sample_id, coverage_sub$sample_id),
    ,
    drop = FALSE
  ]
  if (anyNA(coverage_sub$sample_id)) {
    stop("Missing coverage rows for ", author_label_col, ".")
  }

  if (author_label_col == coverage_label_cols[1]) {
    section_summary$geo_sample_id <- coverage_sub$geo_sample_id
    section_summary$n_spots <- coverage_sub$n_spots
  } else if (any(coverage_sub$n_spots != section_summary$n_spots)) {
    stop("Coverage spot counts differ among label sources.")
  }
  if (any(coverage_sub$author_sample_id != section_summary$author_sample_id)) {
    stop("Author sample IDs do not match the section metadata.")
  }

  section_summary[[paste0("n_", author_label_col, "_labeled_spots")]] <-
    coverage_sub$n_labeled_spots
  section_summary[[paste0("prop_", author_label_col, "_labeled_spots")]] <-
    coverage_sub$prop_labeled_spots
  section_summary[[paste0("n_", author_label_col, "_labels")]] <-
    coverage_sub$n_author_labels
}

biology_sub <- sample_biology[
  match(section_summary$sample_id, sample_biology$sample_id),
  ,
  drop = FALSE
]
if (anyNA(biology_sub$sample_id)) {
  stop("Missing biological summary rows.")
}
section_summary$fraction_structured <- biology_sub$fraction_structured
section_summary$fraction_mes_hyp <- biology_sub$fraction_mes_hyp

for (author_label_col in bd_label_cols) {
  bd_sub <- author_bd[
    author_bd$author_label_col == author_label_col,
    ,
    drop = FALSE
  ]
  bd_sub <- bd_sub[
    match(section_summary$sample_id, bd_sub$sample_id),
    ,
    drop = FALSE
  ]
  if (anyNA(bd_sub$sample_id) || any(bd_sub$boundary_k != boundary_k)) {
    stop("Missing or incompatible BD rows for ", author_label_col, ".")
  }

  for (bd_col in c(
    "n_labels", "mean_discordant", "boundary_density",
    "expected_boundary_density", "relative_boundary_density"
  )) {
    section_summary[[paste0(author_label_col, "_", bd_col)]] <-
      bd_sub[[bd_col]]
  }
}

section_summary <- section_summary[, c(
  "sample_id", "geo_sample_id", "author_sample_id", "tumor_group",
  "region_group", "section_display", "tumor_section_order",
  "tumor_n_sections", "n_spots",
  "n_org1_labeled_spots", "prop_org1_labeled_spots", "n_org1_labels",
  "n_mp_labeled_spots", "prop_mp_labeled_spots", "n_mp_labels",
  "n_layer_labeled_spots", "prop_layer_labeled_spots", "n_layer_labels",
  "fraction_structured", "fraction_mes_hyp",
  "mp_n_labels", "mp_mean_discordant", "mp_boundary_density",
  "mp_expected_boundary_density", "mp_relative_boundary_density",
  "layer_n_labels", "layer_mean_discordant", "layer_boundary_density",
  "layer_expected_boundary_density", "layer_relative_boundary_density"
)]


################################################################################
# Dominant-metaprogram composition
################################################################################

analysis_ix <- as.character(colData(spe)$sample_id) %in% analysis_sample_ids
metaprogram_levels <- sort(unique(as.character(
  colData(spe)$author_mp[analysis_ix & !is.na(colData(spe)$author_mp)]
)))
if (length(metaprogram_levels) < 2L) {
  stop("Fewer than two dominant metaprogram labels were found.")
}

composition <- do.call(rbind, lapply(seq_len(nrow(section_summary)), function(i) {
  sample_ix <- as.character(colData(spe)$sample_id) ==
    section_summary$sample_id[i]
  labels <- as.character(colData(spe)$author_mp[sample_ix])
  n_mp_labeled_spots <- sum(!is.na(labels))
  if (n_mp_labeled_spots == 0L) {
    stop("No dominant-metaprogram labels for ", section_summary$sample_id[i], ".")
  }
  counts <- as.integer(table(factor(
    labels[!is.na(labels)],
    levels = metaprogram_levels
  )))

  data.frame(
    sample_id = section_summary$sample_id[i],
    geo_sample_id = section_summary$geo_sample_id[i],
    author_sample_id = section_summary$author_sample_id[i],
    tumor_group = section_summary$tumor_group[i],
    region_group = section_summary$region_group[i],
    section_display = section_summary$section_display[i],
    tumor_section_order = section_summary$tumor_section_order[i],
    tumor_n_sections = section_summary$tumor_n_sections[i],
    metaprogram = metaprogram_levels,
    n_spots = section_summary$n_spots[i],
    n_mp_labeled_spots = n_mp_labeled_spots,
    prop_mp_labeled_spots = n_mp_labeled_spots / section_summary$n_spots[i],
    n_metaprogram_spots = counts,
    fraction_metaprogram = counts / n_mp_labeled_spots,
    stringsAsFactors = FALSE
  )
}))

composition_counts <- aggregate(
  n_metaprogram_spots ~ sample_id,
  composition,
  sum
)
composition_fractions <- aggregate(
  fraction_metaprogram ~ sample_id,
  composition,
  sum
)
if (any(
  composition_counts$n_metaprogram_spots[
    match(section_summary$sample_id, composition_counts$sample_id)
  ] != section_summary$n_mp_labeled_spots
) || any(abs(
  composition_fractions$fraction_metaprogram[
    match(section_summary$sample_id, composition_fractions$sample_id)
  ] - 1
) > validation_tolerance)) {
  stop("Dominant-metaprogram compositions do not match their denominators.")
}


################################################################################
# Descriptive within-tumor summaries
################################################################################

tumor_metrics <- data.frame(
  metric = c(
    "fraction_structured",
    "fraction_mes_hyp",
    "mp_relative_boundary_density",
    "layer_relative_boundary_density"
  ),
  metric_display = c(
    "Structured-zone fraction",
    "MES-Hyp dominant-metaprogram fraction",
    "Dominant-metaprogram relative boundary density",
    "State-layer relative boundary density"
  ),
  stringsAsFactors = FALSE
)

tumor_summary <- do.call(rbind, lapply(analysis_tumor_groups, function(tumor) {
  tumor_rows <- section_summary[section_summary$tumor_group == tumor, ]
  do.call(rbind, lapply(seq_len(nrow(tumor_metrics)), function(i) {
    values <- tumor_rows[[tumor_metrics$metric[i]]]
    values <- values[is.finite(values)]
    data.frame(
      tumor_group = tumor,
      tumor_n_sections = nrow(tumor_rows),
      metric = tumor_metrics$metric[i],
      metric_display = tumor_metrics$metric_display[i],
      n_nonmissing_sections = length(values),
      mean_value = mean(values),
      median_value = median(values),
      min_value = min(values),
      max_value = max(values),
      range_value = diff(range(values)),
      stringsAsFactors = FALSE
    )
  }))
}))


################################################################################
# Graph sensitivity analysis
################################################################################

# Primary: construct the 6-NN graph from labeled spots. Sensitivity: construct
# the 6-NN graph from all section spots and retain edges with two labeled
# endpoints. Labeled sources without retained neighbors are reported in the QC
# table and excluded from the sensitivity boundary-density mean.

graph_definitions <- c(
  labeled_only_6nn = "6-NN graph constructed from labeled spots only",
  all_section_6nn_labeled_edges = paste0(
    "6-NN graph constructed from all section spots, then restricted to ",
    "directed edges with two labeled endpoints"
  )
)

qc_rows <- list()
qc_i <- 0L

for (i in seq_len(nrow(section_summary))) {
  sample_id <- section_summary$sample_id[i]
  message("Processing graph sensitivity for sample: ", sample_id)
  sample_ix <- as.character(colData(spe)$sample_id) == sample_id
  coords <- spatialCoords(spe)[sample_ix, , drop = FALSE]
  if (ncol(coords) != 2L || any(!is.finite(coords))) {
    stop("Invalid two-dimensional spatial coordinates for ", sample_id, ".")
  }
  all_section_edges <- make_directed_knn_edges(coords, boundary_k)

  for (author_label_col in bd_label_cols) {
    labels <- as.character(
      colData(spe)[[paste0("author_", author_label_col)]][sample_ix]
    )
    labeled_ix <- !is.na(labels)
    if (sum(labeled_ix) <= boundary_k) {
      stop("Too few labeled spots for ", sample_id, " / ", author_label_col, ".")
    }

    labeled_only_summary <- summarize_labeled_edges(
      make_directed_knn_edges(coords[labeled_ix, , drop = FALSE], boundary_k),
      labels[labeled_ix]
    )
    all_section_summary <- summarize_labeled_edges(all_section_edges, labels)

    saved_relative_bd <- section_summary[[
      paste0(author_label_col, "_relative_boundary_density")
    ]][i]
    if (!isTRUE(all.equal(
      labeled_only_summary$relative_boundary_density,
      saved_relative_bd,
      tolerance = validation_tolerance
    ))) {
      stop(
        "Labeled-only relative BD does not match the canonical result for ",
        sample_id, " / ", author_label_col, "."
      )
    }

    qc_i <- qc_i + 1L
    qc_rows[[qc_i]] <- make_qc_row(
      section_summary[i, ],
      author_label_col,
      names(graph_definitions)[1],
      sum(labeled_ix),
      labeled_only_summary
    )
    qc_i <- qc_i + 1L
    qc_rows[[qc_i]] <- make_qc_row(
      section_summary[i, ],
      author_label_col,
      names(graph_definitions)[2],
      nrow(coords),
      all_section_summary
    )
  }
}

qc_summary <- do.call(rbind, qc_rows)


################################################################################
# Save outputs
################################################################################

parameters <- data.frame(
  dataset_id = dataset_id,
  input_files = paste(basename(input_files), collapse = ";"),
  sections_analyzed = nrow(section_metadata),
  tumors_analyzed = length(analysis_tumor_groups),
  tumor_groups = paste(analysis_tumor_groups, collapse = ","),
  boundary_k = boundary_k,
  bd_label_cols = paste(bd_label_cols, collapse = ","),
  composition_label_col = "mp",
  composition_denominator = "spots_with_nonmissing_author_mp",
  primary_graph = graph_definitions[[1]],
  sensitivity_graph = graph_definitions[[2]],
  graph_edges = "directed",
  edge_distance_units = "full-resolution image pixels",
  zero_neighbor_handling =
    "exclude_source_from_bd_mean_and_report_in_qc",
  inferential_tests = "none",
  stringsAsFactors = FALSE
)

write.csv(section_summary, file = section_summary_file, row.names = FALSE)
write.csv(composition, file = composition_file, row.names = FALSE)
write.csv(tumor_summary, file = tumor_summary_file, row.names = FALSE)
write.csv(qc_summary, file = qc_file, row.names = FALSE)
write.csv(parameters, file = parameters_file, row.names = FALSE)
