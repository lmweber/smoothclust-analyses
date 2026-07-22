################################################################################
# Plot Kuppe et al. myocardial infarction boundary density results
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(ggplot2)
library(patchwork)

dataset_id <- "Kuppe"
results_dir <- file.path("results", dataset_id)
summary_dir <- file.path("plots", dataset_id, "summary_plots")
spatial_dir <- file.path("plots", dataset_id, "spatial_maps")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

bd <- read.csv(
  file.path(results_dir, "Kuppe_celltype_niche_BD_section_summary.csv"),
  stringsAsFactors = FALSE
)
spots <- readRDS(file.path(
  results_dir, "Kuppe_celltype_niche_BD_spot_data.rds"
))

group_levels <- c(
  "myogenic-enriched", "ischaemic-enriched", "fibrotic-enriched"
)
group_labels <- c(
  "Myogenic-enriched", "Ischaemic-enriched", "Fibrotic-enriched"
)
group_palette <- c("#0072B2", "#D55E00", "#009E73")
region_levels <- c(
  "Control", "Remote zone", "Border zone", "Ischaemic zone",
  "Fibrotic zone"
)
region_palette <- c("#0072B2", "#56B4E9", "#E69F00", "#D55E00", "#009E73")

bd$sample_group_key <- bd$sample_group
bd$sample_group <- factor(
  bd$sample_group, levels = group_levels, labels = group_labels
)
bd$region_label <- factor(bd$region_label, levels = region_levels)

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "grey92", color = "grey55")
    )
)

save_plot <- function(plot, filename, width, height) {
  ggsave(paste0(filename, ".png"), plot, width = width, height = height,
         dpi = 300)
  ggsave(paste0(filename, ".pdf"), plot, width = width, height = height)
}

bd_long <- rbind(
  data.frame(
    bd[, c("section_id", "sample_id", "patient_id", "sample_group",
           "region_label")],
    metric = "Raw boundary density", value = bd$boundary_density
  ),
  data.frame(
    bd[, c("section_id", "sample_id", "patient_id", "sample_group",
           "region_label")],
    metric = "Relative boundary density",
    value = bd$relative_boundary_density
  )
)
bd_long$metric <- factor(
  bd_long$metric, c("Raw boundary density", "Relative boundary density")
)
jitter <- position_jitter(width = 0.09, height = 0, seed = 123)

p_group <- ggplot(
  bd_long,
  aes(x = sample_group, y = value, color = sample_group, fill = sample_group)
) +
  geom_boxplot(width = 0.58, alpha = 0.15, outlier.shape = NA) +
  geom_point(
    aes(shape = region_label), position = jitter, size = 2.5, alpha = 0.9
  ) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_color_manual(values = group_palette) +
  scale_fill_manual(values = group_palette) +
  labs(
    title = "Myocardial cell-type-niche boundary density by sample group",
    subtitle = paste(
      "Exact author spot annotations; each point is one Visium section"
    ),
    x = NULL, y = NULL, shape = "Anatomical region"
  ) +
  guides(
    color = "none", fill = "none", shape = guide_legend(nrow = 1)
  ) +
  theme(legend.position = "bottom")

save_plot(
  p_group,
  file.path(summary_dir, "Kuppe_celltype_niche_BD_by_sample_group"),
  9.5, 5.5
)

p_region <- ggplot(
  bd_long,
  aes(x = region_label, y = value, color = region_label, fill = region_label)
) +
  geom_boxplot(width = 0.58, alpha = 0.15, outlier.shape = NA) +
  geom_point(position = jitter, size = 2.5, alpha = 0.9) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_color_manual(values = region_palette) +
  scale_fill_manual(values = region_palette) +
  labs(
    title = "Myocardial cell-type-niche boundary density by tissue region",
    subtitle = "Authors' major histomorphological region labels",
    x = NULL, y = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

save_plot(
  p_region,
  file.path(summary_dir, "Kuppe_celltype_niche_BD_by_tissue_region"),
  10, 5.5
)

# Show the section nearest the median relative BD within each author group.
examples <- do.call(rbind, lapply(split(bd, bd$sample_group_key), function(x) {
  x[which.min(abs(
    x$relative_boundary_density - median(x$relative_boundary_density)
  )), ]
}))
examples <- examples[match(group_levels, examples$sample_group_key), ]

niche_levels <- paste0("ctniche_", 1:9)
niche_palette <- setNames(
  grDevices::hcl.colors(length(niche_levels), palette = "Dark 3"),
  niche_levels
)
spatial_theme <- theme_void(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    plot.subtitle = element_text(hjust = 0.5, size = 9)
  )

spatial_plots <- vector("list", 2 * nrow(examples))
for (i in seq_len(nrow(examples))) {
  result <- examples[i, ]
  spot_sub <- spots[spots$section_id == result$section_id, ]
  spot_sub$author_celltype_niche <- factor(
    spot_sub$author_celltype_niche, levels = niche_levels
  )

  niche_plot <- ggplot(
    spot_sub, aes(x = x, y = y, color = author_celltype_niche)
  ) +
    geom_point(size = 1.15) +
    scale_color_manual(
      values = niche_palette, drop = FALSE, labels = paste("Niche", 1:9)
    ) +
    scale_y_reverse() +
    coord_fixed() +
    labs(
      title = paste0(result$sample_id, " (", result$section_id, ")"),
      subtitle = paste(result$region_label, "|", result$sample_group),
      color = "Cell-type niche"
    ) +
    spatial_theme

  local_plot <- ggplot(
    spot_sub, aes(x = x, y = y, color = local_boundary_density)
  ) +
    geom_point(size = 1.15) +
    scale_color_gradientn(
      colors = c("grey95", "#FDB863", "#B2182B"),
      limits = c(0, 1), na.value = "grey85"
    ) +
    scale_y_reverse() +
    coord_fixed() +
    labs(
      title = "Local boundary density",
      subtitle = sprintf("Relative BD = %.3f", result$relative_boundary_density),
      color = "Local BD"
    ) +
    spatial_theme

  if (i > 1L) {
    niche_plot <- niche_plot + theme(legend.position = "none")
    local_plot <- local_plot + theme(legend.position = "none")
  }
  spatial_plots[[2 * i - 1L]] <- niche_plot
  spatial_plots[[2 * i]] <- local_plot
}

p_spatial <- wrap_plots(spatial_plots, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Representative Kuppe myocardial infarction sections",
    subtitle = "One median-relative-BD section from each author sample group"
  )

save_plot(
  p_spatial,
  file.path(spatial_dir, "Kuppe_celltype_niche_example_spatial_maps"),
  11, 11
)
