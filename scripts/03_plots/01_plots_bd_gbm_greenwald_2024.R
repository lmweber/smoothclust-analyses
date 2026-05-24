################################################################################
# Script to plot boundary density metric results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from root directory of smoothclust-analyses repo

library(ggplot2)

dataset_id <- "gbm_greenwald_2024"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id, "boundary_density")
plots_dir <- file.path(repo_dir, "plots", dataset_id, "boundary_density")

summary_file <- file.path(results_dir, "gbm_boundary_density_summary.csv")


# load data

bd <- read.csv(summary_file)

bd$sample_id <- factor(bd$sample_id, levels = sort(unique(bd$sample_id)))
bd$n_clusters <- as.integer(bd$n_clusters)
bd$label_col <- factor(bd$label_col, levels = paste0("smoothclust_k", sort(unique(bd$n_clusters))))


# common plot settings

theme_set(theme_bw(base_size = 11))

sample_palette <- setNames(
  grDevices::hcl.colors(length(levels(bd$sample_id)), palette = "Dark 3"),
  levels(bd$sample_id)
)


# boundary density by number of clusters

ggplot(bd, aes(x = n_clusters, y = boundary_density, group = sample_id, color = sample_id)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_point(size = 1.4) +
  scale_color_manual(values = sample_palette) +
  labs(x = "Number of clusters", y = "Boundary density", color = "Sample") +
  theme(legend.position = "right")

ggsave(
  file.path(plots_dir, "gbm_bd_by_nclusters.png"),
  width = 8,
  height = 5,
  dpi = 300
)


# mean discordant neighbors by number of clusters

ggplot(bd, aes(x = n_clusters, y = mean_discordant, group = sample_id, color = sample_id)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_point(size = 1.4) +
  scale_color_manual(values = sample_palette) +
  labs(x = "Number of clusters", y = "Mean discordant neighbors", color = "Sample") +
  theme(legend.position = "right")

ggsave(
  file.path(plots_dir, "gbm_mean_discordant_by_nclusters.png"),
  width = 8,
  height = 5,
  dpi = 300
)


# relative boundary density by number of clusters

ggplot(bd, aes(x = n_clusters, y = relative_boundary_density, group = sample_id, color = sample_id)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 0.4, alpha = 0.8, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  scale_color_manual(values = sample_palette) +
  labs(x = "Number of clusters", y = "Relative boundary density", color = "Sample") +
  theme(legend.position = "right")

ggsave(
  file.path(plots_dir, "gbm_relative_bd_by_nclusters.png"),
  width = 8,
  height = 5,
  dpi = 300
)


# heatmap of raw boundary density

ggplot(bd, aes(x = factor(n_clusters), y = sample_id, fill = boundary_density)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", name = "BD") +
  labs(x = "Number of clusters", y = "Sample") +
  theme(panel.grid = element_blank())

ggsave(
  file.path(plots_dir, "gbm_bd_heatmap.png"),
  width = 6,
  height = 6,
  dpi = 300
)


# heatmap of analytic relative boundary density

ggplot(bd, aes(x = factor(n_clusters), y = sample_id, fill = relative_boundary_density)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "viridis", name = "Relative BD", na.value = "gray85") +
  labs(x = "Number of clusters", y = "Sample") +
  theme(panel.grid = element_blank())

ggsave(
  file.path(plots_dir, "gbm_relative_bd_heatmap.png"),
  width = 6,
  height = 6,
  dpi = 300
)


# raw boundary density vs analytic expected boundary density

ggplot(bd, aes(x = expected_boundary_density, y = boundary_density, color = factor(n_clusters))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_brewer(palette = "Dark2", name = "Clusters") +
  coord_equal() +
  labs(x = "Expected boundary density", y = "Observed boundary density") +
  theme(legend.position = "right")

ggsave(
  file.path(plots_dir, "gbm_bd_vs_expected_bd.png"),
  width = 5.5,
  height = 5,
  dpi = 300
)
