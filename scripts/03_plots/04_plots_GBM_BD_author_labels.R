################################################################################
# Script to plot author-derived label boundary density results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from root directory of smoothclust-analyses repo

library(ggplot2)
library(SpatialExperiment)

dataset_id <- "GBM"
map_sample_id <- "zh1019inf"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id, "boundary_density")
plots_dir <- file.path(repo_dir, "plots", dataset_id, "boundary_density")
author_plots_dir <- file.path(plots_dir, "author_labels")
supplement_dir <- file.path(author_plots_dir, "supplementary")

dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)

coverage_file <- file.path(results_dir, "GBM_BD_author_labels_coverage.csv")
author_summary_file <- file.path(results_dir, "GBM_BD_author_labels_summary.csv")
spe_plot_file <- file.path(results_dir, "spe_GBM_author_labels_BD_plot.rds")


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

set_label_factors <- function(x) {
  x$author_label_col <- factor(x$author_label_col, levels = label_levels)
  x$author_label_display <- factor(
    x$author_label_display,
    levels = label_display_levels
  )
  x
}


# load data

coverage <- read.csv(coverage_file)
author_bd <- read.csv(author_summary_file)
spe <- readRDS(spe_plot_file)

label_levels <- c("org1", "org2", "mp")
label_display_levels <- c(
  "org1 zones (3 classes)",
  "org2 zones (5 classes)",
  "Leiden metaprogram labels"
)

coverage <- set_label_factors(coverage)
author_bd <- set_label_factors(author_bd)

author_sample_groups <- data.frame(
  sample_id = c(
    "mgh258",
    "zh1007inf", "zh1007nec",
    "zh1019inf", "zh1019t1",
    "zh8811a", "zh8811b", "zh8812", "zh881inf", "zh881t1",
    "zh916bulk", "zh916inf", "zh916t1"
  ),
  tumor_group = c(
    "MGH258",
    "ZH1007", "ZH1007",
    "ZH1019", "ZH1019",
    "ZH881", "ZH881", "ZH881", "ZH881", "ZH881",
    "ZH916", "ZH916", "ZH916"
  ),
  region_group = c(
    "Tumor, unqualified",
    "Infiltrating", "Necrotic",
    "Infiltrating", "T1 contrast-enhancing",
    "Bulk tumor", "Bulk tumor", "Bulk tumor", "Infiltrating", "T1 contrast-enhancing",
    "Bulk tumor", "Infiltrating", "T1 contrast-enhancing"
  ),
  stringsAsFactors = FALSE
)

author_bd <- merge(author_bd, author_sample_groups, by = "sample_id", sort = FALSE)

sample_order <- aggregate(boundary_density ~ sample_id, author_bd, mean)
sample_order <- sample_order[order(sample_order$boundary_density, sample_order$sample_id), ]
sample_levels <- sample_order$sample_id

coverage_order <- aggregate(prop_labeled_spots ~ sample_id, coverage, max)
coverage_order <- coverage_order[order(coverage_order$prop_labeled_spots, coverage_order$sample_id), ]
coverage$sample_id <- factor(coverage$sample_id, levels = coverage_order$sample_id)
author_bd$sample_id <- factor(author_bd$sample_id, levels = sample_levels)

author_bd$tumor_group <- factor(
  author_bd$tumor_group,
  levels = c("MGH258", "ZH1007", "ZH1019", "ZH881", "ZH916")
)
author_bd$region_group <- factor(
  author_bd$region_group,
  levels = c(
    "Bulk tumor",
    "Infiltrating",
    "Necrotic",
    "T1 contrast-enhancing",
    "Tumor, unqualified"
  )
)

# plot settings

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
)

tumor_palette <- setNames(
  grDevices::hcl.colors(length(levels(author_bd$tumor_group)), palette = "Dark 3"),
  levels(author_bd$tumor_group)
)

region_palette <- setNames(
  grDevices::hcl.colors(length(levels(author_bd$region_group)), palette = "Set 2"),
  levels(author_bd$region_group)
)

sample_shapes <- setNames(
  c(16, 17, 15, 18, 0, 1, 2, 3, 4, 5, 6, 8, 10)[seq_along(sample_levels)],
  sample_levels
)


################################################################################
# Supplementary plots: raw boundary density
################################################################################

# Representative author-derived label maps

map_ix <- as.character(colData(spe)$sample_id) == map_sample_id
map_coords <- as.data.frame(spatialCoords(spe)[map_ix, , drop = FALSE])

map_df <- rbind(
  data.frame(
    x = map_coords[, 1],
    y = map_coords[, 2],
    author_label_display = "org1 zones (3 classes)",
    label = colData(spe)$author_org1[map_ix],
    stringsAsFactors = FALSE
  ),
  data.frame(
    x = map_coords[, 1],
    y = map_coords[, 2],
    author_label_display = "org2 zones (5 classes)",
    label = colData(spe)$author_org2[map_ix],
    stringsAsFactors = FALSE
  ),
  data.frame(
    x = map_coords[, 1],
    y = map_coords[, 2],
    author_label_display = "Leiden metaprogram labels",
    label = colData(spe)$author_mp[map_ix],
    stringsAsFactors = FALSE
  )
)

map_df <- map_df[!is.na(map_df$label), ]
map_df$author_label_display <- factor(
  map_df$author_label_display,
  levels = label_display_levels
)

map_labels <- sort(unique(map_df$label))
map_palette <- setNames(
  grDevices::hcl.colors(length(map_labels), palette = "Dark 3"),
  map_labels
)

map_plot <- ggplot(map_df, aes(x = x, y = -y, color = label)) +
  geom_point(size = 0.5) +
  coord_fixed() +
  facet_wrap(~ author_label_display, ncol = 3) +
  scale_color_manual(values = map_palette) +
  labs(
    x = NULL,
    y = NULL,
    color = "Label"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "gray85", color = "gray40"),
    strip.text = element_text(size = 10),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(0.35, "lines")
  )

save_plot(
  map_plot,
  file.path(supplement_dir, paste0("GBM_author_label_raw_BD_maps_", map_sample_id, "_supp")),
  width = 9,
  height = 4
)


# Raw boundary density by label source and sample

raw_bd_sample_plot <- ggplot(
  author_bd,
  aes(x = sample_id, y = boundary_density)
) +
  geom_point(size = 2.2, color = "gray25") +
  facet_wrap(~ author_label_display, nrow = 1, scales = "free_y") +
  labs(
    x = "Sample",
    y = "Raw boundary density"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "none"
  )

save_plot(
  raw_bd_sample_plot,
  file.path(supplement_dir, "GBM_author_label_raw_BD_by_sample_supp"),
  width = 10,
  height = 3.8
)


# Raw boundary density by surgical region

raw_bd_region_plot <- ggplot(
  author_bd,
  aes(
    x = region_group,
    y = boundary_density,
    color = tumor_group,
    shape = sample_id
  )
) +
  stat_summary(
    aes(x = region_group, y = boundary_density, group = region_group),
    inherit.aes = FALSE,
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.35,
    color = "gray35",
    show.legend = FALSE
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 1001),
    size = 2.1,
    alpha = 0.9
  ) +
  facet_wrap(~ author_label_display, nrow = 1, scales = "free_y") +
  scale_color_manual(values = tumor_palette) +
  scale_shape_manual(values = sample_shapes, na.translate = FALSE) +
  labs(
    x = "Surgical region",
    y = "Raw boundary density",
    color = "Tumor",
    shape = "Sample"
  ) +
  guides(
    color = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 3, byrow = TRUE)
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 7)
  )

save_plot(
  raw_bd_region_plot,
  file.path(supplement_dir, "GBM_author_label_raw_BD_by_region_group_supp"),
  width = 10.5,
  height = 5.2
)


# Raw boundary density by tumor group

raw_bd_tumor_plot <- ggplot(
  author_bd,
  aes(x = tumor_group, y = boundary_density, color = region_group)
) +
  stat_summary(
    aes(group = tumor_group),
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.35,
    color = "gray35",
    show.legend = FALSE
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 1001),
    size = 2.1,
    alpha = 0.9
  ) +
  facet_wrap(~ author_label_display, ncol = 1, scales = "free_y") +
  scale_color_manual(values = region_palette) +
  labs(
    x = "Tumor",
    y = "Raw boundary density",
    color = "Surgical region"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )

save_plot(
  raw_bd_tumor_plot,
  file.path(supplement_dir, "GBM_author_label_raw_BD_by_tumor_group_supp"),
  width = 7,
  height = 9
)


################################################################################
# Supplementary plots: relative BD and supporting summaries
################################################################################

# Relative boundary density by sample

relative_bd_sample_plot <- ggplot(
  author_bd,
  aes(x = relative_boundary_density, y = sample_id)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_point(size = 2.2, color = "gray25") +
  facet_wrap(~ author_label_display, ncol = 1, scales = "free_x") +
  labs(
    x = "Relative boundary density",
    y = "Sample"
  ) +
  theme(legend.position = "none")

save_plot(
  relative_bd_sample_plot,
  file.path(supplement_dir, "GBM_author_label_relative_BD_by_sample_supp"),
  width = 7,
  height = 8
)


# Relative boundary density by surgical region

relative_bd_region_plot <- ggplot(
  author_bd,
  aes(x = region_group, y = relative_boundary_density, color = tumor_group)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  stat_summary(
    aes(group = region_group),
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.35,
    color = "gray35",
    show.legend = FALSE
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 1001),
    size = 2.1,
    alpha = 0.9
  ) +
  facet_wrap(~ author_label_display, ncol = 1, scales = "free_y") +
  scale_color_manual(values = tumor_palette) +
  labs(
    x = "Surgical region",
    y = "Relative boundary density",
    color = "Tumor"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )

save_plot(
  relative_bd_region_plot,
  file.path(supplement_dir, "GBM_author_label_relative_BD_by_region_group_supp"),
  width = 7,
  height = 9
)


# Relative boundary density by tumor group

relative_bd_tumor_plot <- ggplot(
  author_bd,
  aes(x = tumor_group, y = relative_boundary_density, color = region_group)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  stat_summary(
    aes(group = tumor_group),
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.35,
    color = "gray35",
    show.legend = FALSE
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 1001),
    size = 2.1,
    alpha = 0.9
  ) +
  facet_wrap(~ author_label_display, ncol = 1, scales = "free_y") +
  scale_color_manual(values = region_palette) +
  labs(
    x = "Tumor",
    y = "Relative boundary density",
    color = "Surgical region"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right"
  )

save_plot(
  relative_bd_tumor_plot,
  file.path(supplement_dir, "GBM_author_label_relative_BD_by_tumor_group_supp"),
  width = 7,
  height = 9
)


# Author label coverage

coverage$has_author_labels <- coverage$n_labeled_spots > 0

coverage_plot <- ggplot(
  coverage,
  aes(x = sample_id, y = prop_labeled_spots, fill = has_author_labels)
) +
  geom_col(width = 0.75) +
  coord_flip() +
  facet_wrap(~ author_label_display, ncol = 1) +
  scale_fill_manual(values = c("FALSE" = "gray75", "TRUE" = "#1B9E77"), guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "Sample",
    y = "Spots with author labels"
  )

save_plot(
  coverage_plot,
  file.path(supplement_dir, "GBM_author_label_coverage_supp"),
  width = 6.5,
  height = 8
)
