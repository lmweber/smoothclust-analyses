################################################################################
# Script to plot smoothing parameter sensitivity boundary density results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(ggplot2)
library(ggspavis)
library(patchwork)
library(SpatialExperiment)

dataset_id <- "GBM"
map_sample_id <- "zh1019t1"
map_n_clusters <- 6
map_smoothing_method <- "uniform"
dose_response_n_clusters <- 6

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)
plots_dir <- file.path(repo_dir, "plots", dataset_id)
sensitivity_plots_dir <- file.path(plots_dir, "sensitivity_smoothing_parameter")
summary_plots_dir <- file.path(sensitivity_plots_dir, "summary_plots")
spatial_maps_dir <- file.path(sensitivity_plots_dir, "spatial_maps")
figure_dir <- file.path(sensitivity_plots_dir, "figure")

dir.create(summary_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_maps_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

bd_summary_file <- file.path(
  results_dir,
  "GBM_BD_sensitivity_smoothing_parameter_summary.csv"
)
rank_stability_file <- file.path(
  results_dir,
  "GBM_BD_sensitivity_smoothing_parameter_rank_stability.csv"
)
spe_plot_file <- file.path(
  results_dir,
  "spe_GBM_BD_sensitivity_smoothing_parameter_plot.rds"
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

parameter_label <- function(parameter_type, value) {
  ifelse(
    parameter_type == "bandwidth",
    sprintf("%.3f", value),
    as.character(as.integer(value))
  )
}

parameter_labeller <- as_labeller(c(
  uniform = 'atop("uniform", "Bandwidth " * italic(b))',
  kernel = 'atop("kernel", "Bandwidth " * italic(b))',
  knn = 'atop("kNN", "Smoothing neighbors " * italic(k)[s])'
), default = label_parsed)

make_spatial_plot <- function(spe_map, map_row) {
  label_col <- map_row$label_col
  n_labels <- map_row$n_labels
  parameter_text <- if (map_row$smoothing_parameter_type == "bandwidth") {
    paste0("b=", sprintf("%.3f", map_row$smoothing_parameter_value))
  } else {
    paste0("k_s=", as.integer(map_row$smoothing_parameter_value))
  }
  
  plotCoords(
    spe_map,
    annotate = label_col,
    pal = cluster_palette[seq_len(n_labels)],
    point_size = 1.15
  ) +
    labs(
      title = parameter_text,
      subtitle = sprintf("Raw BD = %.3f", map_row$boundary_density),
      color = "Cluster"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}


# load data

bd <- read.csv(bd_summary_file)
rank_stability <- read.csv(rank_stability_file)
spe <- readRDS(spe_plot_file)

method_levels <- c("uniform", "kernel", "knn")
method_display <- c(
  uniform = "uniform",
  kernel = "kernel",
  knn = "knn"
)

bd$smoothing_method <- factor(
  bd$smoothing_method,
  levels = method_levels,
  labels = method_display[method_levels]
)
bd$n_clusters <- as.integer(bd$n_clusters)
bd$sample_id <- factor(bd$sample_id, levels = sort(unique(bd$sample_id)))
bd$smoothing_parameter_label <- parameter_label(
  bd$smoothing_parameter_type,
  bd$smoothing_parameter_value
)

rank_stability$smoothing_method <- factor(
  rank_stability$smoothing_method,
  levels = method_levels,
  labels = method_display[method_levels]
)
rank_stability$n_clusters <- as.integer(rank_stability$n_clusters)
rank_stability$smoothing_parameter_label <- parameter_label(
  rank_stability$smoothing_parameter_type,
  rank_stability$smoothing_parameter_value
)

bandwidth_levels <- sprintf(
  "%.3f",
  sort(unique(bd$smoothing_parameter_value[bd$smoothing_parameter_type == "bandwidth"]))
)
knn_levels <- as.character(sort(unique(
  as.integer(bd$smoothing_parameter_value[bd$smoothing_parameter_type == "smoothing_k"])
)))
parameter_levels <- c(bandwidth_levels, knn_levels)
bd$smoothing_parameter_label <- factor(
  bd$smoothing_parameter_label,
  levels = parameter_levels
)
rank_stability$smoothing_parameter_label <- factor(
  rank_stability$smoothing_parameter_label,
  levels = parameter_levels
)
rank_stability$n_clusters_factor <- factor(
  rank_stability$n_clusters,
  levels = sort(unique(rank_stability$n_clusters))
)


# common plot settings

theme_set(
  theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
)

cluster_palette <- c(
  "#0072B2",
  "#E69F00",
  "#009E73",
  "#CC79A7",
  "#56B4E9",
  "#F0E442",
  "#D55E00",
  "#000000"
)


################################################################################
# Illustrative uniform smoothing maps
################################################################################

spe_map <- spe[, as.character(colData(spe)$sample_id) == map_sample_id]

map_bd <- bd[
  bd$sample_id == map_sample_id &
    bd$smoothing_method == map_smoothing_method &
    bd$n_clusters == map_n_clusters,
]

map_small <- map_bd[which.min(map_bd$smoothing_parameter_value), ]
map_large <- map_bd[which.max(map_bd$smoothing_parameter_value), ]

small_bandwidth_map <- make_spatial_plot(spe_map, map_small)
large_bandwidth_map <- make_spatial_plot(spe_map, map_large)

spatial_maps_plot <- small_bandwidth_map + large_bandwidth_map +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

save_plot(
  spatial_maps_plot,
  file.path(
    spatial_maps_dir,
    paste0(
      "GBM_sensitivity_smoothing_parameter_uniform_maps_",
      map_sample_id,
      "_k",
      map_n_clusters
    )
  ),
  width = 8.4,
  height = 4.5
)


################################################################################
# Smoothing dose-response curves
################################################################################

dose_bd <- bd[bd$n_clusters == dose_response_n_clusters, ]
default_lines <- unique(
  dose_bd[
    dose_bd$is_default,
    c("smoothing_method", "smoothing_parameter_value")
  ]
)

dose_response_plot <- ggplot(
  dose_bd,
  aes(
    x = smoothing_parameter_value,
    y = boundary_density,
    group = sample_id
  )
) +
  geom_line(color = "gray65", linewidth = 0.25, alpha = 0.45) +
  geom_point(color = "gray45", size = 0.9, alpha = 0.55) +
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "line",
    linewidth = 0.8,
    color = "black"
  ) +
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "point",
    size = 1.8,
    color = "black"
  ) +
  geom_vline(
    data = default_lines,
    aes(xintercept = smoothing_parameter_value),
    linetype = "dashed",
    color = "gray35",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ smoothing_method, nrow = 1, scales = "free_x",
    labeller = parameter_labeller
  ) +
  labs(
    x = "Smoothing parameter value",
    y = "Raw boundary density",
    subtitle = paste0("K = ", dose_response_n_clusters)
  ) +
  theme(
    strip.text = element_text(size = 11.5),
    plot.subtitle = element_text(size = 11.5)
  )

save_plot(
  dose_response_plot,
  file.path(
    summary_plots_dir,
    paste0("GBM_sensitivity_smoothing_parameter_dose_response_k", dose_response_n_clusters)
  ),
  width = 9,
  height = 3.8
)


################################################################################
# Sample-rank stability across smoothing parameters
################################################################################

rank_stability_plot <- ggplot(
  rank_stability,
  aes(
    x = smoothing_parameter_label,
    y = n_clusters_factor,
    fill = spearman_rho
  )
) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_tile(
    data = rank_stability[rank_stability$is_default, ],
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ smoothing_method, nrow = 1, scales = "free_x",
    labeller = parameter_labeller
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 1),
    na.value = "gray85",
    name = expression(atop("Spearman", rho))
  ) +
  labs(
    x = "Smoothing parameter value",
    y = "Number of clusters"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(size = 11.5)
  )

save_plot(
  rank_stability_plot,
  file.path(summary_plots_dir, "GBM_sensitivity_smoothing_parameter_rank_stability"),
  width = 8.8,
  height = 3.9
)


################################################################################
# Combined supplementary figure
################################################################################

smoothing_parameter_figure <- wrap_elements(spatial_maps_plot) /
  dose_response_plot /
  rank_stability_plot +
  plot_layout(heights = c(1.35, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

save_plot(
  smoothing_parameter_figure,
  file.path(figure_dir, "GBM_supp_sensitivity_smoothing_parameter"),
  width = 10.5,
  height = 10.5
)
