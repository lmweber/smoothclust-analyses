################################################################################
# Script to plot external spatial clustering boundary density results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from root directory of smoothclust-analyses repo

library(ggplot2)
library(ggspavis)
library(patchwork)
library(SpatialExperiment)

dataset_id <- "GBM"
map_sample_id <- "zh1019t1"
map_n_clusters <- 5
boundary_k <- 6

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id, "boundary_density")
plots_dir <- file.path(repo_dir, "plots", dataset_id, "boundary_density")
external_plots_dir <- file.path(plots_dir, "external_clustering")
summary_plots_dir <- file.path(external_plots_dir, "summary_plots")
cluster_plots_dir <- file.path(external_plots_dir, "cluster_maps")
local_bd_plots_dir <- file.path(external_plots_dir, "local_BD_maps")
figure_2_dir <- file.path(external_plots_dir, "figure_2")

dir.create(summary_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cluster_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(local_bd_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_2_dir, recursive = TRUE, showWarnings = FALSE)

smoothclust_summary_file <- file.path(results_dir, "GBM_BD_summary.csv")
baseline_summary_file <- file.path(results_dir, "GBM_BD_baseline_kmeans_summary.csv")
external_summary_file <- file.path(results_dir, "GBM_BD_external_clustering_summary.csv")
external_spe_plot_file <- file.path(results_dir, "spe_GBM_external_clustering_BD_plot.rds")


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

add_standard_columns <- function(x, method, method_display) {
  if (!"n_clusters_requested" %in% names(x)) {
    x$n_clusters_requested <- x$n_clusters
  }
  if (!"n_labels" %in% names(x)) {
    x$n_labels <- x$n_clusters_requested
  }
  if (!"local_boundary_density_col" %in% names(x)) {
    x$local_boundary_density_col <- paste0("local_bd_", x$label_col)
  }
  
  x$method <- method
  x$method_display <- method_display
  x
}

# load data

bd_smooth <- read.csv(smoothclust_summary_file)
bd_baseline <- read.csv(baseline_summary_file)
bd_external <- read.csv(external_summary_file)
spe_external <- readRDS(external_spe_plot_file)

bd_smooth <- add_standard_columns(
  bd_smooth,
  method = "smoothclust_kmeans",
  method_display = "smoothclust"
)
bd_baseline <- add_standard_columns(
  bd_baseline,
  method = "baseline_unsmoothed_kmeans",
  method_display = "Unsmoothed k-means"
)

external_method_labels <- c(
  bayesspace = "BayesSpace",
  banksy = "Banksy",
  clustsignal = "clustSIGNAL"
)
missing_methods <- setdiff(names(external_method_labels), unique(bd_external$method))
if (length(missing_methods) > 0L) {
  stop("Missing external method results: ", paste(missing_methods, collapse = ", "))
}
bd_external$method_display <- external_method_labels[bd_external$method]

common_cols <- Reduce(
  intersect,
  list(names(bd_smooth), names(bd_baseline), names(bd_external))
)

bd <- rbind(
  bd_smooth[, common_cols],
  bd_baseline[, common_cols],
  bd_external[, common_cols]
)

method_order <- c(
  "smoothclust",
  "Unsmoothed k-means",
  "BayesSpace",
  "Banksy",
  "clustSIGNAL"
)

bd$method_display <- factor(bd$method_display, levels = method_order)
bd$sample_id <- factor(bd$sample_id, levels = sort(unique(bd$sample_id)))
bd$n_clusters_requested <- as.integer(bd$n_clusters_requested)
bd$n_labels <- as.integer(bd$n_labels)

bd_external$method_display <- factor(
  bd_external$method_display,
  levels = method_order
)
bd_external$sample_id <- factor(
  bd_external$sample_id,
  levels = levels(bd$sample_id)
)
bd_external$n_clusters_requested <- as.integer(bd_external$n_clusters_requested)
bd_external$n_labels <- as.integer(bd_external$n_labels)
bd$clustering_setting <- ifelse(
  is.na(bd$n_clusters_requested),
  "default",
  paste0("K=", bd$n_clusters_requested)
)
bd_external$clustering_setting <- ifelse(
  is.na(bd_external$n_clusters_requested),
  "default",
  paste0("K=", bd_external$n_clusters_requested)
)
bd_k_requested <- bd[!is.na(bd$n_clusters_requested), ]


# common plot settings

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
)

method_palette <- c(
  "smoothclust" = "#1B9E77",
  "Unsmoothed k-means" = "#D95F02",
  "BayesSpace" = "#7570B3",
  "Banksy" = "#E7298A",
  "clustSIGNAL" = "#66A61E"
)


################################################################################
# Summary plots
################################################################################

# Raw boundary density by method and requested number of clusters

bd_by_nclusters <- ggplot(
  bd_k_requested,
  aes(
    x = n_clusters_requested,
    y = boundary_density,
    group = interaction(sample_id, method_display),
    color = method_display
  )
) +
  geom_line(linewidth = 0.25, alpha = 0.25) +
  geom_point(size = 1.1, alpha = 0.55) +
  stat_summary(
    aes(group = method_display),
    fun = mean,
    geom = "line",
    linewidth = 1.0
  ) +
  stat_summary(
    aes(group = method_display),
    fun = mean,
    geom = "point",
    size = 2.1
  ) +
  scale_color_manual(values = method_palette, drop = FALSE) +
  labs(
    x = "Requested number of clusters",
    y = "Raw boundary density",
    color = "Method"
  ) +
  theme(legend.position = "right")

save_plot(
  bd_by_nclusters,
  file.path(summary_plots_dir, "GBM_external_clustering_BD_by_nclusters"),
  width = 7.5,
  height = 4.8
)


# Relative boundary density by method and requested number of clusters

relative_bd_by_nclusters <- ggplot(
  bd_k_requested,
  aes(
    x = n_clusters_requested,
    y = relative_boundary_density,
    group = interaction(sample_id, method_display),
    color = method_display
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 0.25, alpha = 0.25, na.rm = TRUE) +
  geom_point(size = 1.1, alpha = 0.55, na.rm = TRUE) +
  stat_summary(
    aes(group = method_display),
    fun = mean,
    geom = "line",
    linewidth = 1.0,
    na.rm = TRUE
  ) +
  stat_summary(
    aes(group = method_display),
    fun = mean,
    geom = "point",
    size = 2.1,
    na.rm = TRUE
  ) +
  scale_color_manual(values = method_palette, drop = FALSE) +
  labs(
    x = "Requested number of clusters",
    y = "Relative boundary density",
    color = "Method"
  ) +
  theme(legend.position = "right")

save_plot(
  relative_bd_by_nclusters,
  file.path(summary_plots_dir, "GBM_external_clustering_relative_BD_by_nclusters"),
  width = 7.5,
  height = 4.8
)


# K = 5 plus default ClustSIGNAL method comparison for candidate Figure 2 panel

bd_k_map <- bd[
  bd$n_clusters_requested == map_n_clusters |
    bd$clustering_setting == "default",
]
method_axis_labels <- vapply(
  levels(bd$method_display),
  function(method) {
    method_settings <- unique(
      as.character(bd_k_map$clustering_setting[bd_k_map$method_display == method])
    )
    method_settings <- method_settings[!is.na(method_settings)]
    method_settings <- sub("^K=", "k = ", method_settings)
    
    paste0(method, "\n(", paste(method_settings, collapse = ", "), ")")
  },
  character(1)
)

bd_method_plot <- ggplot(
  bd_k_map,
  aes(x = method_display, y = boundary_density, color = method_display)
) +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "gray35") +
  geom_jitter(width = 0.12, height = 0, size = 1.5, alpha = 0.75) +
  scale_color_manual(values = method_palette, drop = FALSE) +
  scale_x_discrete(labels = method_axis_labels) +
  labs(
    x = NULL,
    y = "Raw boundary density"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, lineheight = 0.9)
  )

save_plot(
  bd_method_plot,
  file.path(figure_2_dir, "GBM_Figure2_external_clustering_raw_BD_by_method"),
  width = 6,
  height = 4.2
)


# Observed label counts for clustSIGNAL

label_count_plot <- ggplot(
  bd_external,
  aes(
    x = clustering_setting,
    y = n_labels,
    color = method_display
  )
) +
  geom_point(size = 1.5, alpha = 0.75) +
  facet_wrap(~ method_display, nrow = 1) +
  scale_color_manual(values = method_palette, drop = FALSE) +
  labs(
    x = "Clustering setting",
    y = "Observed number of labels",
    color = "Method"
  ) +
  theme(legend.position = "none")

save_plot(
  label_count_plot,
  file.path(summary_plots_dir, "GBM_external_clustering_observed_label_counts"),
  width = 8,
  height = 3.2
)


################################################################################
# Representative maps for external methods
################################################################################

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
local_bd_palette <- c("gray92", "#e31a1c")
spe_map <- spe_external[
  ,
  as.character(colData(spe_external)$sample_id) == map_sample_id
]

external_methods <- c("bayesspace", "banksy", "clustsignal")
external_map_plots <- list()
external_local_plots <- list()

for (method in external_methods) {
  default_label_col <- paste0(method, "_default")
  label_col <- if (default_label_col %in% names(colData(spe_map))) {
    default_label_col
  } else {
    paste0(method, "_k", map_n_clusters)
  }
  local_col <- paste0("local_bd_", label_col)
  method_display <- external_method_labels[[method]]
  
  selected_bd <- bd_external[
    bd_external$sample_id == map_sample_id &
      bd_external$method == method &
      (
        bd_external$n_clusters_requested == map_n_clusters |
          is.na(bd_external$n_clusters_requested)
      ) &
      bd_external$boundary_k == boundary_k,
  ]
  
  if (!label_col %in% names(colData(spe_map))) {
    stop("Missing external label column: ", label_col)
  }
  if (!local_col %in% names(colData(spe_map))) {
    stop("Missing local boundary density column: ", local_col)
  }
  
  n_labels <- length(unique(colData(spe_map)[[label_col]]))
  label_palette <- cluster_palette[seq_len(n_labels)]
  
  external_map_plots[[method]] <- plotCoords(
    spe_map,
    annotate = label_col,
    pal = label_palette,
    point_size = 1.2
  ) +
    labs(
      title = method_display,
      subtitle = sprintf(
        "%s; labels = %s; BD = %.3f",
        selected_bd$clustering_setting,
        selected_bd$n_labels,
        selected_bd$boundary_density
      ),
      color = "Label"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 7.5, hjust = 0.5),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
  
  external_local_plots[[method]] <- plotCoords(
    spe_map,
    annotate = local_col,
    pal = local_bd_palette,
    point_size = 1.2
  ) +
    labs(
      title = paste(method_display, "local BD"),
      subtitle = sprintf("k = %s spatial neighbors", boundary_k),
      color = "Local BD"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 7.5, hjust = 0.5),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
  
  save_plot(
    external_map_plots[[method]],
    file.path(
      cluster_plots_dir,
      paste0("GBM_external_cluster_map_", label_col, "_", map_sample_id)
    ),
    width = 4,
    height = 4
  )
  
  save_plot(
    external_local_plots[[method]],
    file.path(
      local_bd_plots_dir,
      paste0("GBM_external_local_BD_map_", label_col, "_", map_sample_id)
    ),
    width = 4.4,
    height = 4
  )
}

external_maps_figure <- (
  wrap_plots(external_map_plots, nrow = 1) /
    wrap_plots(external_local_plots, nrow = 1)
) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_plot(
  external_maps_figure,
  file.path(
    figure_2_dir,
    paste0("GBM_Figure2_external_clustering_maps_", map_sample_id)
  ),
  width = 10.5,
  height = 6.8
)


################################################################################
# Supplementary heatmap-style summary
################################################################################

sample_mean_bd <- aggregate(boundary_density ~ sample_id, bd_k_map, mean)
sample_order <- sample_mean_bd$sample_id[order(sample_mean_bd$boundary_density)]
bd_k_map$sample_id <- factor(bd_k_map$sample_id, levels = rev(sample_order))

bd_heatmap <- ggplot(
  bd_k_map,
  aes(x = method_display, y = sample_id, fill = boundary_density)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_viridis_c(option = "magma", name = "BD") +
  labs(
    x = NULL,
    y = "Sample"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

save_plot(
  bd_heatmap,
  file.path(
    summary_plots_dir,
    paste0("GBM_external_clustering_BD_heatmap_k", map_n_clusters)
  ),
  width = 6.4,
  height = 5.8
)
