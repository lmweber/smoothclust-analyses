################################################################################
# Script to plot boundary density metric results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from root directory of smoothclust-analyses repo

library(ggplot2)
library(ggspavis)
library(SpatialExperiment)

dataset_id <- "GBM"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id, "boundary_density")
plots_dir <- file.path(repo_dir, "plots", dataset_id, "boundary_density")
summary_plots_dir <- file.path(plots_dir, "summary_plots")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_plots_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(results_dir, "GBM_BD_summary.csv")
spe_plot <- file.path(results_dir, "spe_GBM_BD_plot.rds")


# load data

bd <- read.csv(summary_file)

bd$sample_id <- factor(bd$sample_id, levels = sort(unique(bd$sample_id)))
bd$n_clusters <- as.integer(bd$n_clusters)
bd$label_col <- factor(bd$label_col, levels = paste0("smoothclust_k", sort(unique(bd$n_clusters))))


# common plot settings

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
)

sample_palette <- setNames(
  grDevices::hcl.colors(length(levels(bd$sample_id)), palette = "Dark 3"),
  levels(bd$sample_id)
)

local_bd_palette <- c("gray92", "#e31a1c")


# boundary density by number of clusters

bd_nclusters_plot <- ggplot(
  bd,
  aes(x = n_clusters, y = boundary_density, group = sample_id, color = sample_id)
) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_point(size = 1.4) +
  scale_color_manual(values = sample_palette) +
  labs(x = "Number of clusters", y = "Boundary density", color = "Sample") +
  theme(legend.position = "right")

ggsave(
  file.path(summary_plots_dir, "GBM_BD_by_nclusters.png"),
  bd_nclusters_plot,
  width = 8,
  height = 5,
  dpi = 300
)


# mean discordant neighbors by number of clusters

mean_discordant_plot <- ggplot(
  bd,
  aes(x = n_clusters, y = mean_discordant, group = sample_id, color = sample_id)
) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_point(size = 1.4) +
  scale_color_manual(values = sample_palette) +
  labs(x = "Number of clusters", y = "Mean discordant neighbors", color = "Sample") +
  theme(legend.position = "right")

ggsave(
  file.path(summary_plots_dir, "GBM_mean_discordant_by_nclusters.png"),
  mean_discordant_plot,
  width = 8,
  height = 5,
  dpi = 300
)


# relative boundary density by number of clusters

relative_bd_nclusters_plot <- ggplot(
  bd,
  aes(
    x = n_clusters,
    y = relative_boundary_density,
    group = sample_id,
    color = sample_id
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 0.4, alpha = 0.8, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  scale_color_manual(values = sample_palette) +
  labs(x = "Number of clusters", y = "Relative boundary density", color = "Sample") +
  theme(legend.position = "right")

ggsave(
  file.path(summary_plots_dir, "GBM_relative_BD_by_nclusters.png"),
  relative_bd_nclusters_plot,
  width = 8,
  height = 5,
  dpi = 300
)


# heatmap of raw boundary density

bd_heatmap_plot <- ggplot(
  bd,
  aes(x = factor(n_clusters), y = sample_id, fill = boundary_density)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", name = "BD") +
  labs(x = "Number of clusters", y = "Sample") +
  theme(panel.grid = element_blank())

ggsave(
  file.path(summary_plots_dir, "GBM_BD_heatmap.png"),
  bd_heatmap_plot,
  width = 6,
  height = 6,
  dpi = 300
)


# heatmap of analytic relative boundary density

relative_bd_heatmap_plot <- ggplot(
  bd,
  aes(x = factor(n_clusters), y = sample_id, fill = relative_boundary_density)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "viridis", name = "Relative BD", na.value = "gray85") +
  labs(x = "Number of clusters", y = "Sample") +
  theme(panel.grid = element_blank())

ggsave(
  file.path(summary_plots_dir, "GBM_relative_BD_heatmap.png"),
  relative_bd_heatmap_plot,
  width = 6,
  height = 6,
  dpi = 300
)


# raw boundary density vs analytic expected boundary density

expected_bd_plot <- ggplot(
  bd,
  aes(
    x = expected_boundary_density,
    y = boundary_density,
    color = factor(n_clusters)
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_brewer(palette = "Dark2", name = "Clusters") +
  coord_equal() +
  labs(x = "Expected boundary density", y = "Observed boundary density") +
  theme(legend.position = "right")

ggsave(
  file.path(summary_plots_dir, "GBM_BD_vs_expected_BD.png"),
  expected_bd_plot,
  width = 5.5,
  height = 5,
  dpi = 300
)


# spatial maps

spe <- readRDS(spe_plot)

cluster_plots_dir <- file.path(plots_dir, "cluster_maps")
local_bd_plots_dir <- file.path(plots_dir, "local_BD_maps")

dir.create(cluster_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(local_bd_plots_dir, recursive = TRUE, showWarnings = FALSE)

sample_ids <- sort(unique(as.character(colData(spe)$sample_id)))
label_cols <- as.character(levels(bd$label_col))
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
  n_clusters <- as.integer(sub("smoothclust_k", "", label_col))
  pal <- cluster_palette[seq_len(n_clusters)]
  
  for (sample_id in sample_ids) {
    spe_sub <- spe[, as.character(colData(spe)$sample_id) == sample_id]
    
    p <- plotCoords(
      spe_sub,
      annotate = label_col,
      pal = pal,
      point_size = 1.2
    ) +
      labs(color = "Cluster") +
      ggtitle(paste(sample_id, label_col))
    
    ggsave(
      file.path(cluster_plots_dir, paste0("GBM_cluster_map_", label_col, "_", sample_id, ".png")),
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
      file.path(local_bd_plots_dir, paste0("GBM_local_BD_map_", label_col, "_", sample_id, ".png")),
      p,
      width = 4,
      height = 4,
      dpi = 300
    )
  }
}
