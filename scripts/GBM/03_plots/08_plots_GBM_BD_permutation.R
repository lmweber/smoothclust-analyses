################################################################################
# Script to plot permutation-adjusted boundary density results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(ggplot2)
library(patchwork)

dataset_id <- "GBM"
representative_sample_id <- "zh1019t1"
representative_n_clusters <- 6L

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)
plots_dir <- file.path(repo_dir, "plots", dataset_id)
permutation_plots_dir <- file.path(plots_dir, "permutation")
summary_plots_dir <- file.path(permutation_plots_dir, "summary_plots")
figure_dir <- file.path(permutation_plots_dir, "figure")

dir.create(summary_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

bd_summary_file <- file.path(results_dir, "GBM_BD_permutation_summary.csv")
analytic_summary_file <- file.path(results_dir, "GBM_BD_summary.csv")
representative_null_file <- file.path(
  results_dir,
  paste0(
    "GBM_BD_permutation_null_",
    representative_sample_id,
    "_k",
    representative_n_clusters,
    ".csv"
  )
)


# helper functions

save_plot <- function(plot, filename, width, height) {
  ggsave(
    paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
  ggsave(
    paste0(filename, ".pdf"),
    plot = plot,
    width = width,
    height = height
  )
}


# load data

bd <- read.csv(bd_summary_file)
analytic_bd <- read.csv(analytic_summary_file)
representative_null <- read.csv(representative_null_file)

bd$n_clusters <- as.integer(bd$n_clusters)
bd$boundary_k <- as.integer(bd$boundary_k)
bd$sample_id <- factor(bd$sample_id, levels = sort(unique(bd$sample_id)))

analytic_bd$n_clusters <- as.integer(analytic_bd$n_clusters)
analytic_bd$boundary_k <- as.integer(analytic_bd$boundary_k)

analytic_bd <- analytic_bd[
  analytic_bd$n_clusters %in% unique(bd$n_clusters),
  c(
    "sample_id",
    "geo_sample_id",
    "label_col",
    "n_clusters",
    "n_spots",
    "boundary_k",
    "expected_boundary_density",
    "relative_boundary_density"
  )
]

bd_analytic <- merge(
  bd,
  analytic_bd,
  by = c(
    "sample_id",
    "geo_sample_id",
    "label_col",
    "n_clusters",
    "n_spots",
    "boundary_k"
  )
)

representative_result <- bd[
  as.character(bd$sample_id) == representative_sample_id &
    bd$n_clusters == representative_n_clusters,
]

if (nrow(representative_result) != 1L) {
  stop("Expected one representative result row")
}

representative_analytic_result <- bd_analytic[
  bd_analytic$sample_id == representative_sample_id &
    bd_analytic$n_clusters == representative_n_clusters,
]

if (nrow(representative_analytic_result) != 1L) {
  stop("Expected one representative analytic comparison row")
}


# common plot settings

theme_set(
  theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
)

observed_color <- "#D81B60"
permutation_color <- "#2C7FB8"
reference_color <- "gray35"

max_bd <- max(
  bd$boundary_density,
  bd$permutation_mean,
  representative_null$permutation_boundary_density,
  na.rm = TRUE
)

max_expected_bd <- max(
  bd_analytic$expected_boundary_density,
  bd_analytic$permutation_mean,
  na.rm = TRUE
)


################################################################################
# Representative permutation null distribution
################################################################################

null_reference_lines <- data.frame(
  statistic = factor(
    c("Observed BD", "Permutation mean"),
    levels = c("Observed BD", "Permutation mean")
  ),
  value = c(
    representative_result$boundary_density,
    representative_result$permutation_mean
  )
)

null_distribution_plot <- ggplot(
  representative_null,
  aes(x = permutation_boundary_density)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 35,
    fill = permutation_color,
    color = "white",
    alpha = 0.75,
    linewidth = 0.15
  ) +
  geom_vline(
    data = null_reference_lines,
    aes(xintercept = value, color = statistic, linetype = statistic),
    linewidth = 0.75
  ) +
  scale_color_manual(
    values = c("Observed BD" = observed_color, "Permutation mean" = reference_color),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c("Observed BD" = "solid", "Permutation mean" = "dashed"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(0, max_bd * 1.03)) +
  labs(
    title = "Representative permutation null",
    subtitle = sprintf(
      "%s, K = %s, %s permutations",
      representative_sample_id,
      representative_n_clusters,
      representative_result$n_permutations
    ),
    x = "Permuted raw boundary density",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    legend.position = "bottom"
  )

save_plot(
  null_distribution_plot,
  file.path(
    summary_plots_dir,
    paste0(
      "GBM_permutation_null_distribution_",
      representative_sample_id,
      "_k",
      representative_n_clusters
    )
  ),
  width = 5.2,
  height = 4.0
)


################################################################################
# Observed raw BD versus permutation-expected BD
################################################################################

observed_expected_plot <- ggplot(
  bd,
  aes(x = permutation_mean, y = boundary_density)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = reference_color
  ) +
  geom_segment(
    aes(
      xend = permutation_mean,
      yend = permutation_mean
    ),
    color = "gray75",
    linewidth = 0.25
  ) +
  geom_point(
    shape = 21,
    size = 2.4,
    stroke = 0.25,
    fill = permutation_color,
    color = "gray20",
    alpha = 0.85
  ) +
  geom_point(
    data = representative_result,
    shape = 21,
    size = 3.2,
    stroke = 0.55,
    fill = observed_color,
    color = "black"
  ) +
  coord_equal(xlim = c(0, max_bd * 1.03), ylim = c(0, max_bd * 1.03)) +
  labs(
    title = "Observed versus permutation-expected BD",
    subtitle = sprintf(
      "GBM/glioma sections, K = %s, %s permutations",
      representative_n_clusters,
      representative_result$n_permutations
    ),
    x = "Permutation mean BD",
    y = "Observed raw BD"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    legend.position = "none"
  )

save_plot(
  observed_expected_plot,
  file.path(
    summary_plots_dir,
    paste0("GBM_observed_vs_permutation_expected_BD_k", representative_n_clusters)
  ),
  width = 4.7,
  height = 4.0
)


################################################################################
# Analytic expected BD versus permutation-expected BD
################################################################################

analytic_permutation_plot <- ggplot(
  bd_analytic,
  aes(x = expected_boundary_density, y = permutation_mean)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = reference_color
  ) +
  geom_point(
    shape = 21,
    size = 2.4,
    stroke = 0.25,
    fill = permutation_color,
    color = "gray20",
    alpha = 0.85
  ) +
  geom_point(
    data = representative_analytic_result,
    shape = 21,
    size = 3.2,
    stroke = 0.55,
    fill = observed_color,
    color = "black"
  ) +
  coord_equal(
    xlim = c(0, max_expected_bd * 1.03),
    ylim = c(0, max_expected_bd * 1.03)
  ) +
  labs(
    title = "Analytic versus permutation expectation",
    subtitle = sprintf(
      "GBM/glioma sections, K = %s, %s permutations",
      representative_n_clusters,
      representative_result$n_permutations
    ),
    x = "Analytic expected BD",
    y = "Permutation mean BD"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    legend.position = "none"
  )

save_plot(
  analytic_permutation_plot,
  file.path(
    summary_plots_dir,
    paste0("GBM_analytic_vs_permutation_expected_BD_k", representative_n_clusters)
  ),
  width = 4.7,
  height = 4.0
)


################################################################################
# Combined supplementary figure
################################################################################

permutation_figure <- 
  null_distribution_plot + observed_expected_plot + analytic_permutation_plot +
  plot_layout(ncol = 3, widths = c(1.05, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_plot(
  permutation_figure,
  file.path(figure_dir, "GBM_supp_permutation_BD"),
  width = 14.4,
  height = 4.2
)
