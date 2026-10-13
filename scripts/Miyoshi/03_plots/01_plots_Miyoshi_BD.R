################################################################################
# Plot section-level summaries of the Miyoshi et al. annotations
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(ggplot2)
library(patchwork)

dataset_id <- "Miyoshi"
results_dir <- file.path("results", dataset_id)
summary_dir <- file.path("plots", dataset_id, "summary_plots")
spatial_dir <- file.path("plots", dataset_id, "spatial_maps")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

bd <- read.csv(file.path(
  results_dir, "Miyoshi_BD_section_summary.csv"
), stringsAsFactors = FALSE)
counts <- read.csv(file.path(
  results_dir, "Miyoshi_author_label_counts.csv"
), stringsAsFactors = FALSE)
spots <- readRDS(file.path(results_dir, "Miyoshi_BD_spot_data.rds"))

save_plot <- function(plot, filename, width, height) {
  ggsave(paste0(filename, ".png"), plot, width = width, height = height,
         dpi = 300, bg = "white")
  ggsave(paste0(filename, ".pdf"), plot, width = width, height = height,
         bg = "white")
}
theme_set(theme_bw(base_size = 11))

region_levels <- c(
  "L1", "L2/3", "L3/4", "L3/4/5", "L5/6", "L6b", "WM1", "WM2", "WM3"
)
region_palette <- setNames(
  grDevices::hcl.colors(length(region_levels), "Spectral"), region_levels
)

# Reuse the nine-region relative-BD order across section-level plots.
q9 <- bd[bd$label_scheme == "author_cortical_region", ]
q9 <- q9[order(q9$relative_boundary_density), ]
section_levels <- q9$section_id
q9$section_id <- factor(q9$section_id, section_levels)
composition <- counts[counts$label_scheme == "author_cortical_region", ]
composition$section_id <- factor(composition$section_id, section_levels)
composition$label <- factor(composition$label, region_levels)
composition$fraction <- ave(
  composition$n_spots, composition$section_id, FUN = function(x) x / sum(x)
)

p_rank <- ggplot(q9, aes(section_id, relative_boundary_density)) +
  geom_segment(
    aes(xend = section_id, y = 0, yend = relative_boundary_density),
    color = "grey75", linewidth = 0.35
  ) +
  geom_point(aes(color = relative_boundary_density), size = 2.3) +
  scale_color_viridis_c(option = "C", end = 0.9) +
  labs(
    title = "Published cortical-region annotations summarized by section",
    subtitle = "Sections are ordered by relative boundary density",
    x = NULL, y = "Relative boundary density", color = "Relative BD"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
p_composition <- ggplot(
  composition, aes(section_id, fraction, fill = label)
) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = region_palette, drop = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Label composition", fill = "Published region") +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )
p_ranked <- p_rank / p_composition +
  plot_layout(heights = c(1.35, 1))
save_plot(
  p_ranked, file.path(summary_dir, "Miyoshi_BD_ranked_sections"),
  12, 8
)

paired <- bd[, c(
  "section_id", "label_scheme", "relative_boundary_density"
)]
paired$label_scheme <- factor(
  paired$label_scheme,
  c("author_cortical_region", "author_network_region"),
  c("Nine published regions", "WM1-WM3 collapsed")
)
p_paired <- ggplot(
  paired,
  aes(label_scheme, relative_boundary_density, group = section_id)
) +
  geom_line(color = "grey72", linewidth = 0.4) +
  geom_point(aes(color = label_scheme), size = 2.1) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  labs(
    title = "Boundary density depends on the annotation map",
    subtitle = "Lines connect the two summaries for each of 39 sections",
    x = NULL, y = "Relative boundary density"
  ) +
  guides(color = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
save_plot(
  p_paired, file.path(summary_dir, "Miyoshi_BD_q9_vs_WM_collapsed"),
  6.5, 5.5
)

spatial_theme <- theme_void(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    legend.position = "bottom"
  )
# Two section-pairs per row, with ten sections per page in relative-BD order.
atlas_theme <- spatial_theme + theme(
  plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
  plot.subtitle = element_text(hjust = 0.5, size = 11),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 11)
)
page_rows <- split(seq_len(nrow(q9)), ceiling(seq_len(nrow(q9)) / 10))
all_section_pages <- lapply(seq_along(page_rows), function(page_i) {
  rows <- page_rows[[page_i]]
  page_plots <- unlist(lapply(rows, function(i) {
    result <- q9[i, ]
    spot_sub <- spots[
      spots$section_id == result$section_id &
        spots$label_scheme == "author_cortical_region",
    ]
    spot_sub$author_label <- factor(spot_sub$author_label, region_levels)
    list(
      ggplot(spot_sub, aes(x, y, color = author_label)) +
        geom_point(size = 0.65, show.legend = TRUE) +
        scale_color_manual(
          values = region_palette, limits = region_levels, drop = FALSE,
          guide = if (i == rows[1]) guide_legend(nrow = 1) else "none"
        ) +
        scale_y_reverse() + coord_fixed() +
        labs(
          title = as.character(result$section_id),
          subtitle = sprintf(
            "Regions; rBD = %.3f", result$relative_boundary_density
          ),
          color = "Published region"
        ) + atlas_theme,
      ggplot(spot_sub, aes(x, y, color = local_boundary_density)) +
        geom_point(size = 0.65) +
        scale_color_gradientn(
          colors = c("grey95", "#FDB863", "#B2182B"), limits = c(0, 1),
          guide = if (i == rows[1]) guide_colorbar(
            direction = "horizontal", title.position = "top",
            barwidth = grid::unit(4.5, "cm")
          ) else "none"
        ) +
        scale_y_reverse() + coord_fixed() +
        labs(
          title = "Local raw BD", subtitle = "k = 6", color = "Local raw BD"
        ) + atlas_theme
    )
  }), recursive = FALSE)
  wrap_plots(page_plots, ncol = 4, nrow = 5, guides = "collect") +
    plot_annotation(
      title = sprintf(
        "Miyoshi frontal-cortex sections ordered by relative BD (%d-%d of 39)",
        min(rows), max(rows)
      ),
      subtitle = "Read each row left to right: regions, local raw BD, then the next section",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
    ) &
    theme(legend.position = "bottom", legend.box = "vertical")
})
for (i in seq_along(all_section_pages)) {
  save_plot(
    all_section_pages[[i]],
    file.path(
      spatial_dir,
      sprintf("Miyoshi_cortical_region_BD_all_sections_page%02d", i)
    ),
    12, 13
  )
}
pdf(
  file.path(spatial_dir, "Miyoshi_cortical_region_BD_all_sections.pdf"),
  width = 12, height = 13, onefile = TRUE
)
invisible(lapply(all_section_pages, print))
dev.off()

# Compare both label schemes in the median section of the nine-region ordering.
median_id <- as.character(q9$section_id[ceiling(nrow(q9) / 2)])
map_spots <- spots[spots$section_id == median_id, ]
map_spots$label_scheme <- factor(
  map_spots$label_scheme,
  c("author_cortical_region", "author_network_region"),
  c("Nine published regions", "WM1-WM3 collapsed")
)
map_levels <- c(region_levels, "WM")
map_palette <- c(region_palette, WM = "#5B5B5B")
map_spots$author_label <- factor(map_spots$author_label, map_levels)
p_maps <- ggplot(map_spots, aes(x, y, color = author_label)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = map_palette, drop = FALSE) +
  scale_y_reverse() + coord_fixed() +
  facet_wrap(~label_scheme) +
  labs(
    title = paste("Annotation-map comparison:", median_id),
    color = "Published region"
  ) +
  spatial_theme
save_plot(
  p_maps,
  file.path(spatial_dir, "Miyoshi_q9_vs_WM_collapsed_map"),
  10, 5.5
)
