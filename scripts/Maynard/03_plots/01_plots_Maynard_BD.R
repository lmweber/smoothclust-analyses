################################################################################
# Plot Maynard et al. manual-layer boundary density results
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(ggplot2)
library(patchwork)

dataset_id <- "Maynard"
results_dir <- file.path("results", dataset_id)
summary_dir <- file.path("plots", dataset_id, "summary_plots")
spatial_dir <- file.path("plots", dataset_id, "spatial_maps")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

bd <- read.csv(file.path(
  results_dir, "Maynard_BD_section_summary.csv"
), stringsAsFactors = FALSE)
spots <- readRDS(file.path(results_dir, "Maynard_BD_spot_data.rds"))

save_plot <- function(plot, filename, width, height) {
  ggsave(paste0(filename, ".png"), plot, width = width, height = height,
         dpi = 300, bg = "white")
  ggsave(paste0(filename, ".pdf"), plot, width = width, height = height,
         bg = "white")
}
theme_set(theme_bw(base_size = 14))

# Order donor rows with the 0 um pair before the 300 um pair.
bd <- bd[order(bd$patient_id, bd$position_um, bd$adjacent_replicate), ]
bd$section_label <- factor(
  paste0(bd$section_id, "\n", bd$position_um, " \u00b5m"),
  paste0(bd$section_id, "\n", bd$position_um, " \u00b5m")
)
bd$adjacent_replicate <- factor(bd$adjacent_replicate)
patient_palette <- c(
  Br5292 = "#0072B2", Br5595 = "#D55E00", Br8100 = "#009E73"
)

p_summary <- ggplot(
  bd,
  aes(section_label, relative_boundary_density, group = pair_id)
) +
  geom_line(color = "grey55", linewidth = 0.55) +
  geom_point(
    aes(color = patient_id, shape = adjacent_replicate),
    size = 3
  ) +
  scale_color_manual(values = patient_palette) +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    title = "Manual-layer boundary density across 12 DLPFC sections",
    subtitle = paste(
      "Lines connect directly adjacent 10-µm sections;",
      "the two pairs per donor are 300 µm apart"
    ),
    x = "Visium section and position within donor block",
    y = "Relative boundary density",
    color = "Donor", shape = "Adjacent replicate"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
save_plot(
  p_summary, file.path(summary_dir, "Maynard_BD_adjacent_sections"),
  9.5, 5.8
)

layer_levels <- c("L1", "L2", "L3", "L4", "L5", "L6", "WM")
layer_palette <- setNames(
  c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8",
    "#91BFDB", "#4575B4", "#595959"),
  layer_levels
)
spatial_theme <- theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom"
  )
section_ids <- as.character(bd$section_id)
layer_maps <- local_maps <- vector("list", length(section_ids))
for (i in seq_along(section_ids)) {
  section_id <- section_ids[i]
  result <- bd[i, ]
  spot_sub <- spots[spots$section_id == section_id, ]
  spot_sub$author_layer <- factor(spot_sub$author_layer, layer_levels)
  section_subtitle <- sprintf(
    "%s; %s \u00b5m; rep. %s",
    result$patient_id, result$position_um, result$adjacent_replicate
  )

  layer_maps[[i]] <- ggplot(
    spot_sub, aes(x, y, color = author_layer)
  ) +
    geom_point(size = 0.55, show.legend = TRUE) +
    scale_color_manual(
      values = layer_palette, limits = layer_levels, drop = FALSE,
      guide = guide_legend(nrow = 1, override.aes = list(size = 2.5))
    ) +
    scale_y_reverse() + coord_fixed() +
    labs(
      title = paste("Section", section_id),
      subtitle = section_subtitle, color = "Published layer"
    ) + spatial_theme
  # Colors show local raw BD; subtitles report section-level relative BD.
  local_maps[[i]] <- ggplot(
    spot_sub, aes(x, y, color = local_boundary_density)
  ) +
    geom_point(size = 0.55) +
    scale_color_gradientn(
      colors = c("grey95", "#FDB863", "#B2182B"), limits = c(0, 1),
      guide = guide_colorbar(
        direction = "vertical", title.position = "top",
        barheight = grid::unit(4, "cm")
      )
    ) +
    scale_y_reverse() + coord_fixed() +
    labs(
      title = paste("Section", section_id),
      subtitle = sprintf(
        "%s\nrBD = %.3f", section_subtitle, result$relative_boundary_density
      ),
      color = "Local BD"
    ) + spatial_theme
}

map_design <- "A#B#C#D\nE#F#G#H\nI#J#K#L"
map_widths <- c(1, 0.12, 1, 0.12, 1, 0.12, 1)
p_layer_maps <- wrap_plots(
  layer_maps, design = map_design, widths = map_widths, guides = "collect"
) &
  theme(legend.position = "bottom")
p_layer_maps <- p_layer_maps +
  plot_annotation(
    title = "Published manual-layer annotations for all 12 DLPFC sections",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )
# Preserve the paired summary above all annotation maps on the first atlas page.
p_summary_and_layers <- p_summary / wrap_elements(full = p_layer_maps) +
  plot_layout(heights = c(0.8, 2.2))
save_plot(
  p_summary_and_layers,
  file.path(spatial_dir, "Maynard_all_sections_manual_layer_maps"),
  8.5, 10.2
)

p_local_maps <- wrap_plots(
  local_maps, design = map_design, widths = map_widths, guides = "collect"
) &
  theme(legend.position = "right")
p_local_maps <- p_local_maps +
  plot_annotation(
    title = "Local boundary density for all 12 DLPFC sections",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )
save_plot(
  p_local_maps,
  file.path(spatial_dir, "Maynard_all_sections_local_BD_maps"),
  8.5, 6.5
)
