################################################################################
# Script to plot smoothclust vs non-spatial k-means baseline results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from root directory of smoothclust-analyses repo

library(ggplot2)
library(ggspavis)
library(SpatialExperiment)

dataset_id <- "GBM"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)
plots_dir <- file.path(repo_dir, "plots", dataset_id)
comparison_plots_dir <- file.path(plots_dir, "baseline_kmeans_comparison")
summary_plots_dir <- file.path(comparison_plots_dir, "summary_plots")
cluster_plots_dir <- file.path(comparison_plots_dir, "baseline_cluster_maps")
local_bd_plots_dir <- file.path(comparison_plots_dir, "baseline_local_BD_maps")

dir.create(comparison_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cluster_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(local_bd_plots_dir, recursive = TRUE, showWarnings = FALSE)

smoothclust_summary_file <- file.path(results_dir, "GBM_BD_summary.csv")
baseline_summary_file <- file.path(results_dir, "GBM_BD_baseline_kmeans_summary.csv")
baseline_spe_plot_file <- file.path(results_dir, "spe_GBM_baseline_kmeans_BD_plot.rds")


# load data

bd_smooth <- read.csv(smoothclust_summary_file)
bd_baseline <- read.csv(baseline_summary_file)

bd_smooth$method <- "smoothclust_kmeans"
bd_smooth$clustering_method <- "kmeans"
bd_smooth$assay_name <- "counts_smooth"

common_cols <- intersect(names(bd_smooth), names(bd_baseline))
bd <- rbind(
  bd_smooth[, common_cols],
  bd_baseline[, common_cols]
)

method_labels <- c(
  smoothclust_kmeans = "smoothclust + k-means",
  baseline_unsmoothed_kmeans = "Unsmoothed k-means"
)

bd$method <- factor(bd$method, levels = names(method_labels), labels = method_labels)
bd$sample_id <- factor(bd$sample_id, levels = sort(unique(bd$sample_id)))
bd$n_clusters <- as.integer(bd$n_clusters)

bd_wide <- merge(
  bd_smooth,
  bd_baseline,
  by = c("sample_id", "geo_sample_id", "n_clusters", "n_spots", "boundary_k"),
  suffixes = c("_smoothclust", "_baseline")
)

bd_wide$boundary_density_difference <- 
  bd_wide$boundary_density_baseline - bd_wide$boundary_density_smoothclust
bd_wide$relative_boundary_density_difference <- 
  bd_wide$relative_boundary_density_baseline - 
  bd_wide$relative_boundary_density_smoothclust
bd_wide$mean_discordant_difference <- 
  bd_wide$mean_discordant_baseline - bd_wide$mean_discordant_smoothclust
bd_wide$sample_id <- factor(bd_wide$sample_id, levels = sort(unique(bd_wide$sample_id)))


# common plot settings

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
)

method_palette <- c(
  "smoothclust + k-means" = "#1B9E77",
  "Unsmoothed k-means" = "#D95F02"
)

local_bd_palette <- c("gray92", "#e31a1c")


# boundary density by number of clusters

bd_nclusters_plot <- ggplot(bd, aes(
  x = n_clusters,
  y = boundary_density,
  group = interaction(sample_id, method),
  color = method
)) +
  geom_line(linewidth = 0.35, alpha = 0.35) +
  geom_point(size = 1.3, alpha = 0.75) +
  stat_summary(aes(group = method), fun = mean, geom = "line", linewidth = 1.1) +
  stat_summary(aes(group = method), fun = mean, geom = "point", size = 2.2) +
  scale_color_manual(values = method_palette) +
  labs(x = "Number of clusters", y = "Boundary density", color = "Method") +
  theme(legend.position = "right")

ggsave(
  file.path(summary_plots_dir, "GBM_baseline_comparison_BD_by_nclusters.png"),
  bd_nclusters_plot,
  width = 7,
  height = 4.5,
  dpi = 300
)


# relative boundary density by number of clusters

relative_bd_nclusters_plot <- ggplot(bd, aes(
  x = n_clusters,
  y = relative_boundary_density,
  group = interaction(sample_id, method),
  color = method
)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 0.35, alpha = 0.35, na.rm = TRUE) +
  geom_point(size = 1.3, alpha = 0.75, na.rm = TRUE) +
  stat_summary(aes(group = method), fun = mean, geom = "line", linewidth = 1.1, na.rm = TRUE) +
  stat_summary(aes(group = method), fun = mean, geom = "point", size = 2.2, na.rm = TRUE) +
  scale_color_manual(values = method_palette) +
  labs(x = "Number of clusters", y = "Relative boundary density", color = "Method") +
  theme(legend.position = "right")

ggsave(
  file.path(summary_plots_dir, "GBM_baseline_comparison_relative_BD_by_nclusters.png"),
  relative_bd_nclusters_plot,
  width = 7,
  height = 4.5,
  dpi = 300
)


# paired smoothclust vs baseline scatterplots

bd_scatter_plot <- ggplot(
  bd_wide,
  aes(x = boundary_density_smoothclust, y = boundary_density_baseline)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(aes(color = factor(n_clusters)), size = 2, alpha = 0.85) +
  scale_color_brewer(palette = "Dark2", name = "Clusters") +
  coord_equal() +
  labs(
    x = "smoothclust boundary density",
    y = "Unsmoothed baseline boundary density"
  ) +
  theme(legend.position = "right")

ggsave(
  file.path(summary_plots_dir, "GBM_baseline_vs_smoothclust_BD_scatter.png"),
  bd_scatter_plot,
  width = 5.5,
  height = 5,
  dpi = 300
)

relative_bd_scatter_plot <- ggplot(bd_wide, aes(
  x = relative_boundary_density_smoothclust,
  y = relative_boundary_density_baseline
)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(aes(color = factor(n_clusters)), size = 2, alpha = 0.85, na.rm = TRUE) +
  scale_color_brewer(palette = "Dark2", name = "Clusters") +
  coord_equal() +
  labs(
    x = "smoothclust relative boundary density",
    y = "Unsmoothed baseline relative boundary density"
  ) +
  theme(legend.position = "right")

ggsave(
  file.path(summary_plots_dir, "GBM_baseline_vs_smoothclust_relative_BD_scatter.png"),
  relative_bd_scatter_plot,
  width = 5.5,
  height = 5,
  dpi = 300
)


# paired differences

bd_difference_plot <- ggplot(
  bd_wide,
  aes(x = factor(n_clusters), y = boundary_density_difference)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, height = 0, size = 1.3, alpha = 0.7) +
  labs(
    x = "Number of clusters",
    y = "Boundary density difference\n(unsmoothed - smoothclust)"
  )

ggsave(
  file.path(summary_plots_dir, "GBM_baseline_minus_smoothclust_BD_by_nclusters.png"),
  bd_difference_plot,
  width = 6,
  height = 4.5,
  dpi = 300
)

relative_bd_difference_plot <- ggplot(
  bd_wide,
  aes(x = factor(n_clusters), y = relative_boundary_density_difference)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_boxplot(outlier.shape = NA, width = 0.55, na.rm = TRUE) +
  geom_jitter(width = 0.12, height = 0, size = 1.3, alpha = 0.7, na.rm = TRUE) +
  labs(
    x = "Number of clusters",
    y = "Relative boundary density difference\n(unsmoothed - smoothclust)"
  )

ggsave(
  file.path(summary_plots_dir, "GBM_baseline_minus_smoothclust_relative_BD_by_nclusters.png"),
  relative_bd_difference_plot,
  width = 6,
  height = 4.5,
  dpi = 300
)


# heatmap of raw boundary density differences

bd_difference_heatmap <- ggplot(bd_wide, aes(
  x = factor(n_clusters),
  y = sample_id,
  fill = boundary_density_difference
)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#1B9E77",
    mid = "white",
    high = "#D95F02",
    midpoint = 0,
    name = "BD diff."
  ) +
  labs(x = "Number of clusters", y = "Sample") +
  theme(panel.grid = element_blank())

ggsave(
  file.path(summary_plots_dir, "GBM_baseline_minus_smoothclust_BD_heatmap.png"),
  bd_difference_heatmap,
  width = 6,
  height = 6,
  dpi = 300
)


# baseline spatial maps

spe_baseline <- readRDS(baseline_spe_plot_file)

sample_ids <- sort(unique(as.character(colData(spe_baseline)$sample_id)))
label_cols <- paste0("baseline_kmeans_k", sort(unique(bd_baseline$n_clusters)))
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

for (label_col in label_cols) {
  local_col <- paste0("local_bd_", label_col)
  n_clusters <- as.integer(sub("baseline_kmeans_k", "", label_col))
  pal <- cluster_palette[seq_len(n_clusters)]
  
  for (sample_id in sample_ids) {
    spe_sub <- spe_baseline[, as.character(colData(spe_baseline)$sample_id) == sample_id]
    
    p <- plotCoords(
      spe_sub,
      annotate = label_col,
      pal = pal,
      point_size = 1.2
    ) +
      labs(color = "Cluster") +
      ggtitle(paste(sample_id, label_col))
    
    ggsave(
      file.path(cluster_plots_dir, paste0("GBM_baseline_cluster_map_", label_col, "_", sample_id, ".png")),
      p,
      width = 4,
      height = 4,
      dpi = 300
    )
    
    p <- plotCoords(
      spe_sub,
      annotate = local_col,
      pal = local_bd_palette,
      point_size = 1.2
    ) +
      ggtitle(paste(sample_id, local_col))
    
    ggsave(
      file.path(local_bd_plots_dir, paste0("GBM_baseline_local_BD_map_", label_col, "_", sample_id, ".png")),
      p,
      width = 4,
      height = 4,
      dpi = 300
    )
  }
}
