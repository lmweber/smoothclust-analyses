################################################################################
# Script to plot smoothing method sensitivity boundary density results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(ggplot2)
library(patchwork)

dataset_id <- "GBM"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)
plots_dir <- file.path(repo_dir, "plots", dataset_id)
sensitivity_plots_dir <- file.path(plots_dir, "sensitivity_smoothing_method")
summary_plots_dir <- file.path(sensitivity_plots_dir, "summary_plots")
figure_dir <- file.path(sensitivity_plots_dir, "figure")

dir.create(summary_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

bd_summary_file <- file.path(
  results_dir,
  "GBM_BD_sensitivity_smoothing_method_summary.csv"
)
ari_summary_file <- file.path(
  results_dir,
  "GBM_BD_sensitivity_smoothing_method_ari.csv"
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

cluster_labeller <- labeller(
  n_clusters = function(x) paste0("K = ", x)
)


# load data

bd <- read.csv(bd_summary_file)
ari <- read.csv(ari_summary_file)

method_levels <- c("uniform", "kernel", "knn")
method_display <- c(
  uniform = "uniform",
  kernel = "kernel",
  knn = "kNN"
)

bd$smoothing_method <- factor(
  bd$smoothing_method,
  levels = method_levels,
  labels = method_display[method_levels]
)
bd$n_clusters <- as.integer(bd$n_clusters)
bd$sample_id <- factor(bd$sample_id, levels = sort(unique(bd$sample_id)))

ari$method_x <- factor(
  ari$method_x,
  levels = method_levels,
  labels = method_display[method_levels]
)
ari$method_y <- factor(
  ari$method_y,
  levels = method_levels,
  labels = method_display[method_levels]
)
ari$n_clusters <- as.integer(ari$n_clusters)
ari$sample_id <- factor(ari$sample_id, levels = levels(bd$sample_id))


# common plot settings

theme_set(
  theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
)

method_palette <- c(
  uniform = "#0072B2",
  kernel = "#D55E00",
  kNN = "#009E73"
)

k_palette <- setNames(
  c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"),
  as.character(sort(unique(ari$n_clusters)))
)


################################################################################
# Summary plots
################################################################################

# Raw boundary density by smoothing method and number of clusters

raw_bd_plot <- ggplot(
  bd,
  aes(x = smoothing_method, y = boundary_density, fill = smoothing_method)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "gray35", alpha = 0.65) +
  geom_jitter(
    shape = 21,
    width = 0.12,
    height = 0,
    size = 1.1,
    stroke = 0.2,
    color = "gray25",
    alpha = 0.75
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 95,
    size = 5.5,
    color = "black"
  ) +
  facet_wrap(~ n_clusters, nrow = 1, labeller = cluster_labeller) +
  scale_fill_manual(values = method_palette, drop = FALSE) +
  labs(
    x = NULL,
    y = "Raw boundary density",
    fill = "Smoothing method"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

save_plot(
  raw_bd_plot,
  file.path(
    summary_plots_dir,
    "GBM_sensitivity_smoothing_method_raw_BD_by_method_k"
  ),
  width = 10.5,
  height = 4.1
)


# Pairwise adjusted Rand index between smoothing methods

ari$method_pair_display <- factor(
  ari$method_pair,
  levels = c("uniform_vs_kernel", "uniform_vs_knn", "kernel_vs_knn"),
  labels = c("kernel vs uniform", "kNN vs uniform", "kNN vs kernel")
)
ari$n_clusters_factor <- factor(
  ari$n_clusters,
  levels = sort(unique(ari$n_clusters))
)

ari_pairwise_plot <- ggplot(
  ari,
  aes(
    x = n_clusters_factor,
    y = ari,
    color = n_clusters_factor
  )
) +
  geom_point(
    position = position_jitter(width = 0.16, height = 0, seed = 1001),
    size = 1.4,
    alpha = 0.75
  ) +
  stat_summary(
    aes(group = n_clusters_factor),
    fun = median,
    geom = "crossbar",
    width = 0.42,
    linewidth = 0.35,
    color = "gray25",
    show.legend = FALSE
  ) +
  facet_wrap(~ method_pair_display, nrow = 1) +
  scale_color_manual(values = k_palette, name = "K") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Number of clusters",
    y = "Adjusted Rand index"
  ) +
  guides(color = "none") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

save_plot(
  ari_pairwise_plot,
  file.path(
    summary_plots_dir,
    "GBM_sensitivity_smoothing_method_ARI_pairwise_by_k"
  ),
  width = 9,
  height = 4.5
)


# Combined supplementary figure

smoothing_method_figure <- raw_bd_plot / ari_pairwise_plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

save_plot(
  smoothing_method_figure,
  file.path(
    figure_dir,
    "GBM_supp_sensitivity_smoothing_method"
  ),
  width = 10.5,
  height = 7.4
)
