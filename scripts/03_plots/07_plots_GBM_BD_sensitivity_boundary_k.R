################################################################################
# Script to plot BD neighborhood size sensitivity boundary density results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from root directory of smoothclust-analyses repo

library(ggplot2)
library(ggspavis)
library(patchwork)
library(SpatialExperiment)

dataset_id <- "GBM"
map_sample_id <- "zh1019t1"
map_n_clusters <- 6
map_boundary_k_values <- c(2L, 6L, 18L)
default_boundary_k <- 6L

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)
plots_dir <- file.path(repo_dir, "plots", dataset_id)
sensitivity_plots_dir <- file.path(plots_dir, "sensitivity_boundary_k")
summary_plots_dir <- file.path(sensitivity_plots_dir, "summary_plots")
local_maps_dir <- file.path(sensitivity_plots_dir, "local_BD_maps")
figure_dir <- file.path(sensitivity_plots_dir, "figure")

dir.create(summary_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(local_maps_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

bd_summary_file <- file.path(
  results_dir,
  "GBM_BD_sensitivity_boundary_k_summary.csv"
)
spe_plot_file <- file.path(
  results_dir,
  "spe_GBM_BD_sensitivity_boundary_k_plot.rds"
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

boundary_k_label <- function(boundary_k) {
  if (boundary_k == min(map_boundary_k_values)) {
    "Small neighborhood"
  } else if (boundary_k == default_boundary_k) {
    "Default neighborhood"
  } else if (boundary_k == max(map_boundary_k_values)) {
    "Large neighborhood"
  } else {
    paste0("BD k = ", boundary_k)
  }
}

make_local_bd_plot <- function(spe_map, map_row) {
  local_col <- map_row$local_boundary_density_col
  
  plotCoords(
    spe_map,
    annotate = local_col,
    pal = local_bd_palette,
    point_size = 1.15
  ) +
    labs(
      title = boundary_k_label(map_row$boundary_k),
      subtitle = sprintf(
        "BD k = %s; K = %s; BD = %.3f",
        map_row$boundary_k,
        map_row$n_clusters,
        map_row$boundary_density
      ),
      color = "Local BD"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
}


# load data

bd <- read.csv(bd_summary_file)
spe <- readRDS(spe_plot_file)

if (!"label_source" %in% names(bd)) {
  bd$label_source <- "smoothclust"
}

bd <- bd[bd$label_source == "smoothclust", ]
bd$n_clusters <- as.integer(bd$n_clusters)
bd$boundary_k <- as.integer(bd$boundary_k)
bd$sample_id <- factor(bd$sample_id, levels = sort(unique(bd$sample_id)))

boundary_k_levels <- sort(unique(bd$boundary_k))
n_clusters_levels <- sort(unique(bd$n_clusters))

bd$boundary_k_factor <- factor(bd$boundary_k, levels = boundary_k_levels)
bd$n_clusters_factor <- factor(
  bd$n_clusters,
  levels = n_clusters_levels,
  labels = paste0("K = ", n_clusters_levels)
)


# common plot settings

theme_set(
  theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
)

local_bd_palette <- c("gray92", "#e31a1c")
k_palette <- setNames(
  c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"),
  paste0("K = ", n_clusters_levels)
)


################################################################################
# Representative local boundary density maps
################################################################################

spe_map <- spe[, as.character(colData(spe)$sample_id) == map_sample_id]

map_bd <- bd[
  bd$sample_id == map_sample_id &
    bd$n_clusters == map_n_clusters &
    bd$boundary_k %in% map_boundary_k_values,
]
map_bd <- map_bd[match(map_boundary_k_values, map_bd$boundary_k), ]

local_map_plots <- lapply(
  seq_len(nrow(map_bd)),
  function(i) make_local_bd_plot(spe_map, map_bd[i, ])
)

local_maps_plot <- wrap_plots(
  plotlist = local_map_plots,
  nrow = 1,
  guides = "collect"
) &
  theme(legend.position = "right")

save_plot(
  local_maps_plot,
  file.path(
    local_maps_dir,
    paste0(
      "GBM_sensitivity_boundary_k_local_BD_maps_",
      map_sample_id,
      "_k",
      map_n_clusters
    )
  ),
  width = 9.6,
  height = 3.7
)


################################################################################
# Raw BD versus BD neighborhood size
################################################################################

raw_bd_boundary_k_plot <- ggplot(
  bd,
  aes(
    x = boundary_k_factor,
    y = boundary_density,
    fill = n_clusters_factor,
    color = n_clusters_factor
  )
) +
  geom_vline(
    xintercept = which(boundary_k_levels == default_boundary_k),
    linetype = "dashed",
    color = "gray35"
  ) +
  geom_boxplot(
    position = position_dodge2(width = 0.82, preserve = "single"),
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.45,
    linewidth = 0.35
  ) +
  scale_fill_manual(values = k_palette, name = "Clusters") +
  scale_color_manual(values = k_palette, name = "Clusters") +
  labs(
    x = "BD neighborhood size k",
    y = "Raw boundary density"
  ) +
  guides(color = "none") +
  theme(legend.position = "right")

save_plot(
  raw_bd_boundary_k_plot,
  file.path(summary_plots_dir, "GBM_sensitivity_boundary_k_raw_BD_by_boundary_k"),
  width = 7.6,
  height = 4.2
)


################################################################################
# Sample-level raw BD across BD neighborhood sizes
################################################################################

sample_level_bd_plot <- ggplot(
  bd,
  aes(
    x = boundary_k_factor,
    y = boundary_density,
    fill = n_clusters_factor,
    color = n_clusters_factor
  )
) +
  geom_vline(
    xintercept = which(boundary_k_levels == default_boundary_k),
    linetype = "dashed",
    color = "gray35"
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.58,
    alpha = 0.45,
    linewidth = 0.35
  ) +
  geom_jitter(width = 0.12, height = 0, size = 0.55, alpha = 0.35) +
  facet_wrap(~ n_clusters_factor, nrow = 2) +
  scale_fill_manual(values = k_palette, name = "Clusters") +
  scale_color_manual(values = k_palette, name = "Clusters") +
  labs(
    x = "BD neighborhood size k",
    y = "Sample-level raw boundary density"
  ) +
  guides(fill = "none", color = "none") +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

save_plot(
  sample_level_bd_plot,
  file.path(summary_plots_dir, "GBM_sensitivity_boundary_k_sample_level_raw_BD"),
  width = 8.4,
  height = 4.8
)


################################################################################
# Combined supplementary figure
################################################################################

sample_level_bd_panel <- plot_spacer() + sample_level_bd_plot + plot_spacer() +
  plot_layout(widths = c(0.12, 1, 0.12))

boundary_k_figure <- wrap_elements(local_maps_plot) /
  raw_bd_boundary_k_plot /
  sample_level_bd_panel +
  plot_layout(heights = c(1.15, 1, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_plot(
  boundary_k_figure,
  file.path(figure_dir, "GBM_supp_sensitivity_boundary_k"),
  width = 10.5,
  height = 10.6
)
