################################################################################
# Annotation definitions, local raw BD, and paired section summaries
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(ggplot2)
library(patchwork)

results_dir <- file.path("results", "Miyoshi")
summary_dir <- file.path("plots", "Miyoshi", "summary_plots")
bd <- read.csv(file.path(results_dir, "Miyoshi_BD_section_summary.csv"),
               stringsAsFactors = FALSE)
spots <- readRDS(file.path(results_dir, "Miyoshi_BD_spot_data.rds"))
schemes <- c("author_cortical_region", "author_network_region")
tolerance <- 1e-12
region_levels <- c(
  "L1", "L2/3", "L3/4", "L3/4/5", "L5/6", "L6b", "WM1", "WM2", "WM3"
)
region_palette <- setNames(
  grDevices::hcl.colors(length(region_levels), "Spectral"), region_levels
)
map_palette <- c(region_palette, WM = "#5B5B5B")

# Check saved pairs without changing their membership, order, or values.
stopifnot(
  "Summary must contain the required columns" = all(c(
    "section_id", "label_scheme", "n_spots", "n_labeled_spots", "n_labels",
    "boundary_k", "boundary_density", "expected_boundary_density",
    "relative_boundary_density"
  ) %in% names(bd)),
  "Spot data must contain the required columns" = all(c(
    "spot_id", "section_id", "sample_id", "label_scheme", "author_label",
    "x", "y", "local_boundary_density"
  ) %in% names(spots)),
  "Expected 39 complete section pairs" = nrow(bd) == 78L &&
    length(unique(bd$section_id)) == 39L && !anyNA(bd$section_id) &&
    setequal(bd$label_scheme, schemes) &&
    all(table(bd$section_id, bd$label_scheme) == 1L),
  "Spot data must cover the same sections and definitions" =
    !anyNA(spots$section_id) && !anyNA(spots$label_scheme) &&
    setequal(spots$section_id, bd$section_id) &&
    setequal(spots$label_scheme, schemes),
  "All summaries must use k=6 and finite BD values" =
    all(bd$boundary_k == 6L) &&
    all(is.finite(as.matrix(bd[, c(
      "boundary_density", "expected_boundary_density", "relative_boundary_density"
    )]))) && all(bd$expected_boundary_density > 0)
)
for (section_id in unique(bd$section_id)) {
  original <- spots[spots$section_id == section_id &
                      spots$label_scheme == schemes[1], ]
  merged <- spots[spots$section_id == section_id &
                    spots$label_scheme == schemes[2], ]
  section_summary <- bd[bd$section_id == section_id, ]
  raw_bd <- section_summary$boundary_density[
    match(schemes, section_summary$label_scheme)
  ]
  expected_bd <- section_summary$expected_boundary_density[
    match(schemes, section_summary$label_scheme)
  ]
  expected_labels <- ifelse(original$author_label %in% c("WM1", "WM2", "WM3"),
                            "WM", as.character(original$author_label))
  stopifnot(
    "Paired spots must have identical IDs, order, and coordinates" =
      nrow(original) > 6L && !anyNA(original$spot_id) &&
      !anyDuplicated(original$spot_id) &&
      identical(as.character(original$spot_id), as.character(merged$spot_id)) &&
      identical(as.character(original$sample_id), as.character(merged$sample_id)) &&
      identical(unname(as.matrix(original[, c("x", "y")])),
                unname(as.matrix(merged[, c("x", "y")]))),
    "Only WM1-WM3 may be merged" =
      all(original$author_label %in% region_levels) &&
      identical(expected_labels, as.character(merged$author_label)),
    "Merging must not increase local raw BD" =
      all(merged$local_boundary_density <=
            original$local_boundary_density + tolerance),
    "Merging must not increase sample-level raw BD" =
      raw_bd[2] <= raw_bd[1] + tolerance,
    "Merging must not increase expected BD" =
      expected_bd[2] <= expected_bd[1] + tolerance
  )
  for (scheme_i in seq_along(schemes)) {
    local <- if (scheme_i == 1L) original else merged
    saved <- bd[bd$section_id == section_id &
                  bd$label_scheme == schemes[scheme_i], ]
    counts <- as.numeric(table(as.character(local$author_label)))
    n <- nrow(local)
    expected <- 1 - sum(counts * (counts - 1)) / (n * (n - 1))
    stopifnot(
      "Saved spots, counts, and local raw BD must be complete" =
        n == saved$n_spots && n == saved$n_labeled_spots &&
        length(counts) == saved$n_labels &&
        all(is.finite(as.matrix(local[, c("x", "y", "local_boundary_density")]))) &&
        all(local$local_boundary_density >= 0 & local$local_boundary_density <= 1),
      "Saved local, raw, expected, and relative BD must agree" =
        abs(mean(local$local_boundary_density) - saved$boundary_density) <= tolerance &&
        abs(expected - saved$expected_boundary_density) <= tolerance &&
        abs(saved$boundary_density / expected -
              saved$relative_boundary_density) <= tolerance
    )
  }
}
message("Saved BD consistency checks passed for all 39 section pairs.")

# Preserve the existing rule: the middle section in original relative-BD order.
q9 <- bd[bd$label_scheme == schemes[1], ]
q9 <- q9[order(q9$relative_boundary_density), ]
median_rank <- ceiling(nrow(q9) / 2)
median_id <- q9$section_id[median_rank]
stopifnot("Representative section differs from the approved example" =
            median_id == "Nov_24_2021_VisiumHuman_12")
map_spots <- spots[spots$section_id == median_id, ]
x_limits <- range(map_spots$x)
y_limits <- rev(range(map_spots$y))
map_spots$author_label <- factor(map_spots$author_label, names(map_palette))
bd$condition <- factor(bd$label_scheme, schemes, c("Nine regions", "WM merged"))

spatial_theme <- theme_void(base_size = 13) + theme(
  legend.position = "bottom", legend.box = "vertical"
)
spatial_map <- function(scheme, local_bd = FALSE) {
  data <- map_spots[map_spots$label_scheme == scheme, ]
  p <- ggplot(data, aes(x, y)) +
    scale_x_continuous(limits = x_limits, expand = expansion(mult = 0.025)) +
    scale_y_reverse(limits = y_limits, expand = expansion(mult = 0.025)) +
    coord_fixed() + spatial_theme
  if (local_bd) {
    p + geom_point(aes(color = local_boundary_density), size = 0.8) +
      scale_color_gradientn(
        colors = c("grey95", "#FDB863", "#B2182B"), limits = c(0, 1),
        name = NULL, na.value = "grey85",
        guide = guide_colorbar(direction = "vertical",
                               barheight = grid::unit(3, "cm"),
                               barwidth = grid::unit(0.35, "cm"), order = 2)
      ) + theme(legend.position = "right")
  } else {
    p + geom_point(aes(color = author_label), size = 0.8) +
      scale_color_manual(
        values = map_palette, drop = TRUE,
        name = "Region", na.value = "grey85",
        guide = guide_legend(
          nrow = 2, byrow = TRUE, title.position = "top",
          override.aes = list(size = 3), order = 1
        )
      ) + theme(
        legend.title = element_text(hjust = 0.5, margin = margin(t = 8, b = 5)),
        legend.text = element_text(hjust = 0.5)
      )
  }
}
plot_heading <- function(title, size = 4.2, hjust = 0.5, fontface = "plain") {
  ggplot() +
    annotate("text", x = hjust, y = 0, label = title,
             hjust = hjust, fontface = fontface, size = size) +
    xlim(0, 1) + ylim(-1, 1) + theme_void() +
    theme(plot.margin = margin(0, 5.5, 0, 5.5))
}
map_pair <- function(scheme, letter, title) {
  # Align subtitles in one row; keep the BD scale beside its map and Region below.
  wrap_plots(
    A = plot_heading(paste(letter, title, sep = "  "),
                     size = 5.5, hjust = 0, fontface = "bold"),
    B = plot_heading("Labels"), C = plot_heading("Local raw BD"),
    D = spatial_map(scheme),
    E = wrap_plots(spatial_map(scheme, TRUE), guides = "keep"),
    F = guide_area(),
    design = "AA\nBC\nDE\nFF", heights = c(0.4, 0.25, 2.25, 0.8),
    guides = "collect"
  )
}
paired_plot <- function(metric, title, letter) {
  ggplot(bd, aes(condition, .data[[metric]], group = section_id)) +
    geom_line(color = "grey72", linewidth = 0.4) +
    geom_point(aes(color = condition), size = 2.1) +
    scale_color_manual(values = c("#0072B2", "#D55E00"), guide = "none") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = NULL, y = title, tag = letter) + theme_bw(base_size = 14) +
    theme(plot.tag = element_text(face = "bold", size = 16),
          plot.tag.position = "topleft", panel.grid.minor = element_blank())
}

# Top row: existing spatial panels A/B.
# Bottom row: raw BD (C), expected BD (D), and relative BD (E).
# Use separate row layouts so C/D/E have equal widths independent of A/B.
figure <- wrap_plots(
  wrap_plots(
    A = map_pair(schemes[1], "A", "Nine-region annotation"),
    B = map_pair(schemes[2], "B", "WM1-WM3 merged"),
    nrow = 1, widths = c(1, 1), guides = "keep"
  ),
  wrap_plots(
    C = paired_plot("boundary_density", "Raw BD", "C"),
    D = paired_plot("expected_boundary_density", "Expected BD", "D"),
    E = paired_plot("relative_boundary_density", "Relative BD", "E"),
    nrow = 1, widths = c(1, 1, 1), guides = "keep"
  ),
  ncol = 1,
  heights = c(3.4, 2.8),
  guides = "keep"
)

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
filename <- file.path(summary_dir, "Miyoshi_annotation_definition_BD_figure")
ggsave(paste0(filename, ".pdf"), figure, width = 12, height = 7, bg = "white")
ggsave(paste0(filename, ".png"), figure, width = 12, height = 7,
       dpi = 300, bg = "white")
write.csv(data.frame(
  section_id = median_id,
  selection_rule = "Middle section after increasing original-annotation relative BD",
  original_relative_bd_rank = median_rank,
  n_sections = nrow(q9),
  n_retained_spots = q9$n_labeled_spots[median_rank],
  original_relative_boundary_density = q9$relative_boundary_density[median_rank]
), file.path(summary_dir, "Miyoshi_annotation_definition_BD_selection.csv"),
row.names = FALSE)
