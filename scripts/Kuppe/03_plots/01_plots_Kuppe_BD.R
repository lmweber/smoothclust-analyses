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
label_counts <- read.csv(
  file.path(results_dir, "Kuppe_author_celltype_niche_counts.csv"),
  stringsAsFactors = FALSE
)

group_levels <- c(
  "myogenic-enriched", "ischaemic-enriched", "fibrotic-enriched"
)
group_labels <- c(
  "Myogenic-enriched", "Ischaemic-enriched", "Fibrotic-enriched"
)
group_labels_wrapped <- sub("-", "-\n", group_labels, fixed = TRUE)
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
bd$sample_group_index <- match(bd$sample_group_key, group_levels)
bd$x_offset <- ave(
  seq_len(nrow(bd)), bd$sample_group_key,
  FUN = function(i) seq(-0.16, 0.16, length.out = length(i))
)
bd$x_plot <- bd$sample_group_index + bd$x_offset

patient_counts <- table(bd$patient_id)
repeated_patients <- sort(names(patient_counts)[patient_counts > 1])
bd$patient_display <- ifelse(
  bd$patient_id %in% repeated_patients,
  bd$patient_id,
  "Single-section patient"
)
bd$patient_display <- factor(
  bd$patient_display,
  levels = c(repeated_patients, "Single-section patient")
)
patient_shapes <- c(
  setNames(
    c(17, 15, 18, 8, 4)[seq_along(repeated_patients)],
    repeated_patients
  ),
  "Single-section patient" = 16
)

theme_set(
  theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "grey92", color = "grey55")
    )
)

save_plot <- function(plot, filename, width, height) {
  ggsave(paste0(filename, ".png"), plot, width = width, height = height,
         dpi = 300, bg = "white")
  ggsave(paste0(filename, ".pdf"), plot, width = width, height = height,
         bg = "white")
}

jitter <- position_jitter(width = 0.09, height = 0, seed = 123)

p_relative <- ggplot(
  bd,
  aes(x = x_plot, y = relative_boundary_density)
) +
  geom_boxplot(
    aes(
      x = sample_group_index, group = sample_group_key,
      fill = sample_group
    ),
    width = 0.42, alpha = 0.2, outlier.shape = NA, show.legend = FALSE
  ) +
  geom_point(
    aes(shape = patient_display, color = dominant_niche_gt80),
    size = 3.0, alpha = 0.95
  ) +
  scale_x_continuous(
    breaks = seq_along(group_levels), labels = group_labels_wrapped
  ) +
  scale_fill_manual(values = setNames(group_palette, group_labels)) +
  scale_color_manual(
    values = c(`FALSE` = "grey25", `TRUE` = "#E00000"),
    guide = "none"
  ) +
  scale_shape_manual(
    values = patient_shapes, breaks = names(patient_shapes),
    labels = c(repeated_patients, "Single section"), drop = FALSE
  ) +
  labs(
    title = "Relative BD",
    subtitle = "Red: >80% in one niche",
    x = NULL, y = "Relative BD",
    shape = "Patient"
  ) +
  guides(
    shape = guide_legend(
      nrow = 1, title.position = "top",
      override.aes = list(color = "grey25")
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(hjust = 0.5, size = 15),
    legend.text = element_text(size = 13),
    legend.key.width = grid::unit(0.9, "lines"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_expected <- ggplot(
  bd,
  aes(
    x = expected_boundary_density, y = boundary_density,
    color = sample_group
  )
) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = 2) +
  geom_point(
    aes(shape = dominant_niche_gt80), size = 2.8, alpha = 0.9
  ) +
  geom_text(
    data = bd[bd$dominant_niche_gt80, ],
    aes(label = section_id), size = 3.8, color = "black",
    nudge_x = 0.015, nudge_y = 0.025, hjust = 0, show.legend = FALSE
  ) +
  scale_color_manual(
    values = setNames(group_palette, group_labels),
    breaks = group_labels, labels = group_labels_wrapped
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 1),
    labels = c(`FALSE` = "No", `TRUE` = "Yes")
  ) +
  coord_equal() +
  labs(
    title = "Raw vs expected BD",
    subtitle = "Relative BD = raw / expected",
    x = "Expected BD",
    y = "Raw BD",
    color = "Sample category",
    shape = ">80% one niche"
  ) +
  guides(
    color = guide_legend(nrow = 1, title.position = "top"),
    shape = "none"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(hjust = 0.5, size = 15),
    legend.text = element_text(size = 13)
  )

p_group <- p_relative + p_expected +
  plot_layout(widths = c(0.8, 1)) +
  plot_annotation(
    title = "BD of myocardial cell-type niches",
    tag_levels = "A"
  )

save_plot(
  p_group,
  file.path(summary_dir, "Kuppe_celltype_niche_BD_by_sample_group"),
  12.5, 7.5
)

p_region <- ggplot(
  bd,
  aes(
    x = region_label, y = relative_boundary_density,
    color = sample_group, fill = region_label
  )
) +
  geom_boxplot(width = 0.58, alpha = 0.15, outlier.shape = NA) +
  geom_point(position = jitter, size = 2.5, alpha = 0.9) +
  scale_color_manual(values = setNames(group_palette, group_labels)) +
  scale_fill_manual(values = region_palette) +
  labs(
    title = "Relative boundary density by histomorphological region",
    subtitle = paste(
      "Published computational cell-type-niche labels;",
      "each point is one Visium section"
    ),
    x = NULL, y = "Relative boundary density",
    color = "Histomorphological\nsample category"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  guides(fill = "none")

save_plot(
  p_region,
  file.path(summary_dir, "Kuppe_celltype_niche_BD_by_tissue_region"),
  9, 5.5
)

# Display niche composition to make highly imbalanced sections transparent.
label_counts$niche_fraction <- ave(
  label_counts$n_spots, label_counts$section_id,
  FUN = function(x) x / sum(x)
)
composition_metadata <- bd[
  match(label_counts$section_id, bd$section_id),
  c(
    "section_id", "relative_boundary_density", "dominant_niche_gt80"
  )
]
label_counts <- cbind(
  label_counts, composition_metadata[, -1, drop = FALSE]
)
section_order <- bd$section_id[order(
  bd$sample_group_index, bd$relative_boundary_density
)]
section_labels <- sprintf(
  "%s / %s (rBD %.2f)%s",
  bd$section_id, bd$patient_id, bd$relative_boundary_density,
  ifelse(bd$dominant_niche_gt80, " [>80%]", "")
)
names(section_labels) <- bd$section_id
label_counts$section_label <- factor(
  section_labels[label_counts$section_id],
  levels = rev(section_labels[section_order])
)
label_counts$sample_group <- factor(
  label_counts$sample_group, levels = group_levels, labels = group_labels
)

p_composition <- ggplot(
  label_counts,
  aes(
    x = niche_fraction, y = section_label,
    fill = author_celltype_niche
  )
) +
  geom_col(width = 0.82) +
  facet_grid(sample_group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = setNames(
      grDevices::hcl.colors(9, palette = "Dark 3"),
      paste0("ctniche_", 1:9)
    ),
    labels = paste("Niche", 1:9)
  ) +
  scale_x_continuous(labels = scales::percent) +
  labs(
    title = "Composition of published computational cell-type niches",
    subtitle = paste(
      "Sections ordered within sample category by relative BD;",
      "[>80%] denotes >80% one niche"
    ),
    x = "Fraction of labeled spots", y = NULL, fill = "Cell-type niche"
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0)
  )

save_plot(
  p_composition,
  file.path(summary_dir, "Kuppe_celltype_niche_composition_by_section"),
  10, 10
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
spatial_theme <- theme_void(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
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
    geom_point(size = 1.15, show.legend = TRUE) +
    scale_color_manual(
      values = niche_palette, limits = niche_levels, drop = FALSE,
      labels = paste("Niche", 1:9)
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
      title = "Local raw BD",
      subtitle = sprintf("Relative BD = %.3f", result$relative_boundary_density),
      color = "Local raw BD"
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
    title = "Representative Kuppe myocardial infarction sections"
  )

save_plot(
  p_spatial,
  file.path(spatial_dir, "Kuppe_celltype_niche_example_spatial_maps"),
  11, 11
)
