################################################################################
# Script to plot clustering-method recovery of author-derived labels
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(ggplot2)
library(patchwork)
library(scales)
library(SpatialExperiment)

dataset_id <- "GBM"
map_sample_id <- "zh1019inf"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)
plots_dir <- file.path(repo_dir, "plots", dataset_id)
comparison_plots_dir <- file.path(plots_dir, "author_label_clustering_comparison")
supplement_dir <- file.path(comparison_plots_dir, "supplementary")

dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)

comparison_summary_file <- file.path(
  results_dir,
  "GBM_BD_author_label_clustering_comparison_summary.csv"
)
spe_plot_file <- file.path(
  results_dir,
  "spe_GBM_author_label_clustering_comparison_plot.rds"
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

method_k_label <- function(method_display, n_clusters_requested) {
  ifelse(
    is.na(n_clusters_requested),
    paste0(method_display, "\ndefault"),
    paste0(method_display, "\nK = ", n_clusters_requested)
  )
}

label_or_unlabeled <- function(label) {
  label <- as.character(label)
  label[is.na(label)] <- "unlabeled"
  label
}


# load data

comparison_summary <- read.csv(comparison_summary_file)
spe <- readRDS(spe_plot_file)

comparison_summary$target_label_col <- factor(
  comparison_summary$target_label_col,
  levels = c("org1", "org2")
)
comparison_summary$target_label_display <- factor(
  comparison_summary$target_label_display,
  levels = c("org1 zones (3 classes)", "org2 zones (5 classes)")
)
comparison_summary$method <- factor(
  comparison_summary$method,
  levels = c("smoothclust", "bayesspace", "banksy", "clustsignal")
)
comparison_summary$method_display <- factor(
  comparison_summary$method_display,
  levels = c("smoothclust", "BayesSpace", "Banksy", "ClustSIGNAL")
)
comparison_summary$region_group <- factor(
  comparison_summary$region_group,
  levels = c(
    "Bulk tumor",
    "Infiltrating",
    "Necrotic",
    "T1 contrast-enhancing",
    "Tumor, unqualified"
  )
)
comparison_summary$tumor_group <- factor(
  comparison_summary$tumor_group,
  levels = c("MGH258", "ZH1007", "ZH1019", "ZH881", "ZH916")
)

matched_summary <- comparison_summary[comparison_summary$is_matched_granularity, ]
matched_summary$method_k_display <- method_k_label(
  as.character(matched_summary$method_display),
  matched_summary$n_clusters_requested
)

method_k_levels <- unique(matched_summary$method_k_display[
  order(matched_summary$method, matched_summary$n_clusters_requested)
])
matched_summary$method_k_display <- factor(
  matched_summary$method_k_display,
  levels = method_k_levels
)

finite_ari <- matched_summary$adjusted_rand_index[
  is.finite(matched_summary$adjusted_rand_index)
]
if (length(finite_ari) == 0L) {
  stop("No finite ARI values found in matched clustering summary.")
}
ari_y_breaks <- pretty(c(0, finite_ari), n = 5)
ari_y_limits <- range(ari_y_breaks)


# plot settings

theme_set(
  theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "gray90", color = "gray45")
    )
)

method_palette <- setNames(
  grDevices::hcl.colors(length(levels(comparison_summary$method_display)), palette = "Dark 3"),
  levels(comparison_summary$method_display)
)

region_palette <- setNames(
  grDevices::hcl.colors(length(levels(comparison_summary$region_group)), palette = "Set 2"),
  levels(comparison_summary$region_group)
)


################################################################################
# Panel B: ARI by method and target label set
################################################################################

ari_plot <- ggplot(
  matched_summary,
  aes(x = method_k_display, y = adjusted_rand_index, color = region_group)
) +
  stat_summary(
    aes(group = method_k_display),
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.35,
    color = "gray30",
    show.legend = FALSE
  ) +
  geom_point(
    aes(shape = tumor_group),
    position = position_jitter(width = 0.1, height = 0, seed = 1001),
    size = 2.1,
    alpha = 0.9
  ) +
  facet_wrap(~ target_label_display, nrow = 1, scales = "free_x") +
  scale_color_manual(values = region_palette) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(breaks = ari_y_breaks) +
  coord_cartesian(ylim = ari_y_limits) +
  labs(
    x = "Spatial clustering method",
    y = "Adjusted Rand index",
    color = "Surgical region",
    shape = "Tumor"
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  )

save_plot(
  ari_plot,
  file.path(supplement_dir, "GBM_author_label_clustering_ARI_by_method_supp"),
  width = 8.5,
  height = 4.5
)


################################################################################
# Panel C: relationship between spatial compactness and recovery of author labels
################################################################################

bd_ari_plot <- ggplot(
  matched_summary,
  aes(
    x = relative_boundary_density,
    y = adjusted_rand_index,
    color = method_display,
    shape = region_group
  )
) +
  geom_point(size = 2.2, alpha = 0.9) +
  facet_wrap(~ target_label_display, nrow = 1) +
  scale_color_manual(values = method_palette) +
  scale_y_continuous(breaks = ari_y_breaks) +
  coord_cartesian(ylim = ari_y_limits) +
  labs(
    x = "Relative boundary density",
    y = "Adjusted Rand index",
    color = "Method",
    shape = "Surgical region"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

save_plot(
  bd_ari_plot,
  file.path(supplement_dir, "GBM_author_label_clustering_ARI_vs_relative_BD_supp"),
  width = 8.5,
  height = 4.5
)


################################################################################
# Panel A: representative author and method maps
################################################################################

map_summary <- matched_summary[matched_summary$sample_id == map_sample_id, ]

map_ix <- as.character(colData(spe)$sample_id) == map_sample_id
map_coords <- as.data.frame(spatialCoords(spe)[map_ix, , drop = FALSE])
map_df <- data.frame()
panel_levels <- character()

for (target_label_col in levels(matched_summary$target_label_col)) {
  target_summary <- map_summary[
    as.character(map_summary$target_label_col) == target_label_col,
  ]
  if (nrow(target_summary) == 0L) {
    stop("No map summary rows found for target label: ", target_label_col)
  }
  target_display <- unique(as.character(target_summary$target_label_display))
  if (length(target_display) != 1L) {
    stop("Expected one display label for target label: ", target_label_col)
  }
  
  target_labels <- colData(spe)[[paste0("author_", target_label_col)]][map_ix]
  author_panel <- paste0("Author ", target_label_col, "\n", target_display)
  panel_levels <- c(panel_levels, author_panel)
  
  map_df <- rbind(
    map_df,
    data.frame(
      x = map_coords[, 1],
      y = map_coords[, 2],
      panel = author_panel,
      label = label_or_unlabeled(target_labels),
      stringsAsFactors = FALSE
    )
  )
  
  for (method in levels(matched_summary$method)) {
    method_row <- target_summary[
      as.character(target_summary$method) == method,
    ]
    if (nrow(method_row) != 1L) {
      stop(
        "Expected one map row for ",
        target_label_col,
        " and method ",
        method,
        "; found ",
        nrow(method_row)
      )
    }
    
    method_label_col <- as.character(method_row$label_col)
    method_labels <- colData(spe)[[method_label_col]][map_ix]
    method_labels[is.na(target_labels)] <- NA
    method_panel <- paste0(
      as.character(method_row$method_display),
      ifelse(
        is.na(method_row$n_clusters_requested),
        "\ndefault",
        paste0("\nK=", method_row$n_clusters_requested)
      ),
      "; ARI ",
      sprintf("%.2f", method_row$adjusted_rand_index)
    )
    panel_levels <- c(panel_levels, method_panel)
    
    map_df <- rbind(
      map_df,
      data.frame(
        x = map_coords[, 1],
        y = map_coords[, 2],
        panel = method_panel,
        label = label_or_unlabeled(method_labels),
        stringsAsFactors = FALSE
      )
    )
  }
}

map_df$panel <- factor(map_df$panel, levels = panel_levels)

map_labels <- c(sort(setdiff(unique(map_df$label), "unlabeled")), "unlabeled")
map_labels <- map_labels[map_labels %in% unique(map_df$label)]
map_palette <- setNames(
  grDevices::hcl.colors(length(map_labels), palette = "Dark 3"),
  map_labels
)
map_palette["unlabeled"] <- "gray85"

map_plot <- ggplot(map_df, aes(x = x, y = -y, color = label)) +
  geom_point(size = 0.32) +
  coord_fixed() +
  facet_wrap(~ panel, nrow = 2) +
  scale_color_manual(values = map_palette) +
  labs(x = NULL, y = NULL, color = "Label") +
  theme_void(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "gray90", color = "gray45"),
    strip.text = element_text(size = 8),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.key.height = grid::unit(0.3, "lines")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))

save_plot(
  map_plot,
  file.path(supplement_dir, paste0("GBM_author_label_clustering_maps_", map_sample_id, "_supp")),
  width = 13,
  height = 5.8
)


################################################################################
# Supplementary Figure
################################################################################

summary_row_plot <- wrap_plots(
  ari_plot,
  bd_ari_plot,
  nrow = 1
)

s6_plot <- wrap_plots(
  map_plot,
  summary_row_plot,
  ncol = 1,
  heights = c(1.2, 1)
) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")

save_plot(
  s6_plot,
  file.path(supplement_dir, "GBM_author_label_clustering_comparison_supp"),
  width = 13,
  height = 10
)
