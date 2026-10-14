################################################################################
# Script to assemble figure for boundary density metric results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(ggplot2)
library(ggspavis)
library(patchwork)
library(SpatialExperiment)

dataset_id <- "GBM"
sample_id <- "zh1019t1"
n_clusters <- 6
boundary_k <- 6

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)
plots_dir <- file.path(repo_dir, "plots", dataset_id)
comparison_plots_dir <- file.path(plots_dir, "baseline_kmeans_comparison")
figure_dir <- file.path(comparison_plots_dir, "figure")
paper_figures_dir <- file.path(repo_dir, "..", "smoothclust_paper", "figures")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

smoothclust_summary_file <- file.path(results_dir, "GBM_BD_summary.csv")
baseline_summary_file <- file.path(results_dir, "GBM_BD_baseline_kmeans_summary.csv")
smoothclust_spe_plot_file <- file.path(results_dir, "spe_GBM_BD_plot.rds")
baseline_spe_plot_file <- file.path(results_dir, "spe_GBM_baseline_kmeans_BD_plot.rds")

input_files <- c(
  smoothclust_summary_file,
  baseline_summary_file,
  smoothclust_spe_plot_file,
  baseline_spe_plot_file
)


# load data

bd_smooth <- read.csv(smoothclust_summary_file)
bd_baseline <- read.csv(baseline_summary_file)

spe_smooth <- readRDS(smoothclust_spe_plot_file)
spe_baseline <- readRDS(baseline_spe_plot_file)

bd_smooth$method <- "smoothclust_kmeans"
bd_baseline$method <- "baseline_unsmoothed_kmeans"

bd_wide <- merge(
  bd_smooth,
  bd_baseline,
  by = c("sample_id", "geo_sample_id", "n_clusters", "n_spots", "boundary_k"),
  suffixes = c("_smoothclust", "_baseline")
)

bd_wide$boundary_density_difference <- 
  bd_wide$boundary_density_baseline - bd_wide$boundary_density_smoothclust

selected_result <- bd_wide[
  bd_wide$sample_id == sample_id &
    bd_wide$n_clusters == n_clusters &
    bd_wide$boundary_k == boundary_k,
]

bd_smooth$n_clusters <- as.integer(bd_smooth$n_clusters)
bd_baseline$n_clusters <- as.integer(bd_baseline$n_clusters)
bd_wide$n_clusters <- as.integer(bd_wide$n_clusters)

sample_mean_bd <- aggregate(boundary_density ~ sample_id, bd_smooth, mean)
sample_order <- sample_mean_bd$sample_id[order(sample_mean_bd$boundary_density)]

bd_smooth$sample_id <- factor(bd_smooth$sample_id, levels = rev(sample_order))
bd_wide$sample_id <- factor(bd_wide$sample_id, levels = rev(sample_order))


# plot settings

theme_set(
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
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
domain_palette <- cluster_palette[seq_len(n_clusters)]
local_bd_palette <- c("gray92", "#e31a1c")
k_palette <- setNames(
  c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"),
  as.character(3:8)
)

smooth_label_col <- paste0("smoothclust_k", n_clusters)
baseline_label_col <- paste0("baseline_kmeans_k", n_clusters)
local_bd_col <- paste0("local_bd_", smooth_label_col)

spe_smooth_sub <- spe_smooth[, as.character(colData(spe_smooth)$sample_id) == sample_id]
spe_baseline_sub <- spe_baseline[, as.character(colData(spe_baseline)$sample_id) == sample_id]


# spatial domain maps

p_smooth_domains <- plotCoords(
  spe_smooth_sub,
  annotate = smooth_label_col,
  pal = domain_palette,
  point_size = 1.2
) +
  labs(
    title = "smoothclust domains",
    subtitle = sprintf("%s, K = %s, Raw BD = %.3f",
                       sample_id, n_clusters, selected_result$boundary_density_smoothclust),
    color = "Cluster"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5),
    legend.title = element_text(size = 10.5),
    legend.text = element_text(size = 10)
  )

p_local_bd <- plotCoords(
  spe_smooth_sub,
  annotate = local_bd_col,
  pal = local_bd_palette,
  point_size = 1.2
) +
  labs(
    title = "Local raw BD",
    subtitle = sprintf("smoothclust domains, k = %s neighbors", boundary_k),
    color = "Local raw BD"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5),
    legend.title = element_text(size = 10.5),
    legend.text = element_text(size = 10)
  )

p_baseline_domains <- plotCoords(
  spe_baseline_sub,
  annotate = baseline_label_col,
  pal = domain_palette,
  point_size = 1.2
) +
  labs(
    title = "Unsmoothed k-means",
    subtitle = sprintf("%s, K = %s, Raw BD = %.3f",
                       sample_id, n_clusters, selected_result$boundary_density_baseline),
    color = "Cluster"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5),
    legend.title = element_text(size = 10.5),
    legend.text = element_text(size = 10)
  )


# summary panels

selected_point <- bd_wide[
  bd_wide$sample_id == sample_id &
    bd_wide$n_clusters == n_clusters &
    bd_wide$boundary_k == boundary_k,
]

max_bd <- max(
  bd_wide$boundary_density_smoothclust,
  bd_wide$boundary_density_baseline,
  na.rm = TRUE
)

p_scatter <- ggplot(
  bd_wide,
  aes(x = boundary_density_smoothclust, y = boundary_density_baseline)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray45") +
  geom_point(
    aes(fill = factor(n_clusters)), shape = 21,
    size = 1.6, stroke = 0.25, color = "gray45", alpha = 0.75
  ) +
  geom_point(
    data = selected_point,
    shape = 21,
    size = 3.2,
    stroke = 0.7,
    fill = "#D81B60",
    color = "black"
  ) +
  scale_fill_manual(values = k_palette, name = "K") +
  coord_equal(xlim = c(0, max_bd), ylim = c(0, max_bd), expand = TRUE) +
  labs(
    title = "Sample-level raw BD",
    x = "smoothclust raw BD",
    y = "Unsmoothed k-means raw BD"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "right"
  )

p_difference <- ggplot(
  bd_wide,
  aes(x = factor(n_clusters), y = boundary_density_difference)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45") +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = "gray92") +
  geom_jitter(width = 0.12, height = 0, size = 1.1, alpha = 0.55) +
  geom_point(
    data = selected_point,
    shape = 21,
    size = 3.0,
    stroke = 0.7,
    fill = "#D81B60",
    color = "black"
  ) +
  labs(
    title = "Paired raw BD difference",
    x = "Number of clusters",
    y = "Raw BD difference\n(unsmoothed - smoothclust)"
  ) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))

p_heatmap <- ggplot(
  bd_smooth,
  aes(x = factor(n_clusters), y = sample_id, fill = boundary_density)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_tile(
    data = bd_smooth[
      bd_smooth$sample_id == sample_id &
        bd_smooth$n_clusters == n_clusters &
        bd_smooth$boundary_k == boundary_k,
    ],
    fill = NA,
    color = "#D81B60",
    linewidth = 0.9
  ) +
  scale_fill_viridis_c(option = "magma", name = "Raw BD") +
  labs(
    title = "smoothclust raw BD",
    x = "Number of clusters",
    y = "Sample"
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.text.y = element_text(size = 10, margin = margin(r = 4)),
    legend.position = "right"
  )


# assemble and save

fig <- (
  p_smooth_domains + p_local_bd + p_baseline_domains +
    plot_layout(ncol = 3, widths = c(1, 1.15, 1))
) / (
  p_scatter + p_difference + p_heatmap +
    plot_layout(ncol = 3, widths = c(1, 1, 1.2))
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 10.5),
    legend.title = element_text(size = 10.5),
    legend.text = element_text(size = 10)
  )

png_file <- file.path(figure_dir, "GBM_BD_figure_zh1019t1.png")
pdf_file <- file.path(figure_dir, "GBM_BD_figure_zh1019t1.pdf")

ggsave(png_file, fig, width = 10.5, height = 7.2, dpi = 300)
ggsave(pdf_file, fig, width = 10.5, height = 7.2)
