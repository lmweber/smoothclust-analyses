################################################################################
# Script to plot author-derived label boundary density results
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
author_plots_dir <- file.path(plots_dir, "author_labels")
supplement_dir <- file.path(author_plots_dir, "supplementary")

dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)

coverage_file <- file.path(results_dir, "GBM_BD_author_labels_coverage.csv")
author_summary_file <- file.path(results_dir, "GBM_BD_author_labels_summary.csv")
sample_biology_file <- file.path(results_dir, "GBM_BD_author_labels_sample_biology.csv")
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
  x$author_label_short <- factor(x$author_label_short, levels = label_levels)
  x
}

rho_label <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) {
    return("rho = NA")
  }
  rho <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  paste0("rho = ", sprintf("%.2f", rho))
}


# load data

coverage <- read.csv(coverage_file)
author_bd <- read.csv(author_summary_file)
sample_biology <- read.csv(sample_biology_file)
spe <- readRDS(spe_plot_file)

label_levels <- c("org1", "org2", "mp", "layer")
label_display_levels <- c(
  "org1 zones (3 classes)",
  "org2 zones (5 classes)",
  "Dominant metaprogram labels",
  "State layers (5 classes)"
)

expected_mp_display <- "Dominant metaprogram labels"
if (!identical(
  unique(coverage$author_label_display[coverage$author_label_col == "mp"]),
  expected_mp_display
) || !identical(
  unique(author_bd$author_label_display[author_bd$author_label_col == "mp"]),
  expected_mp_display
)) {
  stop(
    "Author-label inputs do not use the expected dominant-metaprogram terminology."
  )
}

coverage <- set_label_factors(coverage)
author_bd <- set_label_factors(author_bd)

author_bd <- merge(
  author_bd,
  sample_biology[, setdiff(names(sample_biology), "n_labeled_spots")],
  by = c("sample_id", "geo_sample_id", "author_sample_id", "n_spots"),
  sort = FALSE
)

sample_order <- aggregate(relative_boundary_density ~ sample_id, author_bd, mean)
sample_order <- sample_order[order(sample_order$relative_boundary_density, sample_order$sample_id), ]
sample_levels <- sample_order$sample_id

coverage_order <- aggregate(prop_labeled_spots ~ sample_id, coverage, max)
coverage_order <- coverage_order[order(coverage_order$prop_labeled_spots, coverage_order$sample_id), ]
coverage$sample_id <- factor(coverage$sample_id, levels = coverage_order$sample_id)
author_bd$sample_id <- factor(author_bd$sample_id, levels = sample_levels)
sample_biology$sample_id <- factor(sample_biology$sample_id, levels = sample_levels)


# plot settings

theme_set(
  theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "gray90", color = "gray45")
    )
)

label_palette <- setNames(
  grDevices::hcl.colors(length(label_levels), palette = "Dark 3"),
  label_levels
)

################################################################################
# Panel A: representative label maps and local BD maps
################################################################################

map_ix <- as.character(colData(spe)$sample_id) == map_sample_id
map_coords <- as.data.frame(spatialCoords(spe)[map_ix, , drop = FALSE])

map_metrics <- author_bd[as.character(author_bd$sample_id) == map_sample_id, ]
map_metrics <- map_metrics[match(label_levels, as.character(map_metrics$author_label_col)), ]
map_metrics$facet_label <- paste0(
  map_metrics$author_label_short,
  "\nBD ", sprintf("%.3f", map_metrics$boundary_density),
  "; relBD ", sprintf("%.2f", map_metrics$relative_boundary_density)
)

label_map_df <- data.frame()
local_map_df <- data.frame()

for (i in seq_len(nrow(map_metrics))) {
  label_col <- as.character(map_metrics$label_col[i])
  local_col <- as.character(map_metrics$local_boundary_density_col[i])
  
  label_map_df <- rbind(
    label_map_df,
    data.frame(
      x = map_coords[, 1],
      y = map_coords[, 2],
      author_label_col = as.character(map_metrics$author_label_col[i]),
      facet_label = map_metrics$facet_label[i],
      label = as.character(colData(spe)[[label_col]][map_ix]),
      stringsAsFactors = FALSE
    )
  )
  
  local_map_df <- rbind(
    local_map_df,
    data.frame(
      x = map_coords[, 1],
      y = map_coords[, 2],
      author_label_col = as.character(map_metrics$author_label_col[i]),
      facet_label = map_metrics$facet_label[i],
      local_boundary_density = colData(spe)[[local_col]][map_ix],
      stringsAsFactors = FALSE
    )
  )
}

label_map_df <- label_map_df[!is.na(label_map_df$label), ]
local_map_df <- local_map_df[!is.na(local_map_df$local_boundary_density), ]

label_map_df$facet_label <- factor(label_map_df$facet_label, levels = map_metrics$facet_label)
local_map_df$facet_label <- factor(local_map_df$facet_label, levels = map_metrics$facet_label)

map_labels <- sort(unique(label_map_df$label))
map_palette <- setNames(
  grDevices::hcl.colors(length(map_labels), palette = "Dark 3"),
  map_labels
)

label_map_plot <- ggplot(label_map_df, aes(x = x, y = -y, color = label)) +
  geom_point(size = 0.35) +
  coord_fixed() +
  facet_wrap(~ facet_label, nrow = 1) +
  scale_color_manual(values = map_palette) +
  labs(x = NULL, y = NULL, color = "Label") +
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "gray90", color = "gray45"),
    strip.text = element_text(size = 10),
    legend.position = "right",
    legend.text = element_text(size = 7.5),
    legend.key.height = grid::unit(0.4, "lines")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))

local_map_plot <- ggplot(
  local_map_df,
  aes(x = x, y = -y, color = local_boundary_density)
) +
  geom_point(size = 0.35) +
  coord_fixed() +
  facet_wrap(~ facet_label, nrow = 1) +
  scale_color_gradientn(
    colors = rev(grDevices::hcl.colors(60, palette = "YlOrRd")),
    limits = c(0, 1)
  ) +
  labs(x = NULL, y = NULL, color = "Local BD") +
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "gray90", color = "gray45"),
    strip.text = element_text(size = 10),
    legend.position = "right"
  )

map_pair_plot <- label_map_plot / local_map_plot +
  plot_layout(heights = c(1, 1), guides = "collect")

save_plot(
  map_pair_plot,
  file.path(supplement_dir, paste0("GBM_author_label_revised_maps_", map_sample_id, "_supp")),
  width = 11,
  height = 5.8
)


################################################################################
# Panel B: sample-level relative BD by author label source
################################################################################

relative_bd_label_plot <- ggplot(
  author_bd,
  aes(x = author_label_display, y = relative_boundary_density, group = sample_id)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray45") +
  geom_line(color = "gray78", linewidth = 0.3) +
  stat_summary(
    aes(group = author_label_display),
    fun = median,
    geom = "crossbar",
    width = 0.55,
    linewidth = 0.35,
    color = "gray25",
    show.legend = FALSE
  ) +
  geom_point(aes(fill = author_label_short), shape = 21, size = 2.2, color = "gray20") +
  scale_fill_manual(values = label_palette, guide = "none") +
  coord_cartesian(ylim = c(0, 1.03)) +
  labs(
    x = "Author label source",
    y = "Relative boundary density"
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  )

save_plot(
  relative_bd_label_plot,
  file.path(supplement_dir, "GBM_author_label_revised_relative_BD_by_label_source_supp"),
  width = 6.5,
  height = 4.2
)


################################################################################
# Panel C: structured organization versus metaprogram/layer BD
################################################################################

biology_label_levels <- c("mp", "layer")
biology_bd <- author_bd[as.character(author_bd$author_label_col) %in% biology_label_levels, ]
biology_bd$author_label_display <- factor(
  biology_bd$author_label_display,
  levels = label_display_levels[label_levels %in% biology_label_levels]
)

rho_df <- do.call(
  rbind,
  lapply(
    split(biology_bd, biology_bd$author_label_display),
    function(x) {
      data.frame(
        author_label_display = unique(x$author_label_display),
        x = 0.03,
        y = 0.96,
        label = rho_label(x$fraction_structured, x$relative_boundary_density),
        stringsAsFactors = FALSE
      )
    }
  )
)

structured_bd_plot <- ggplot(
  biology_bd,
  aes(
    x = fraction_structured,
    y = relative_boundary_density,
    color = fraction_mes_hyp
  )
) +
  geom_point(size = 2.5, alpha = 0.95) +
  geom_text(
    data = rho_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.5,
    color = "gray25"
  ) +
  facet_wrap(~ author_label_display, nrow = 1) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_gradientn(
    colors = grDevices::hcl.colors(60, palette = "BluYl"),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = "Structured-zone fraction",
    y = "Relative boundary density",
    color = "MES-Hyp fraction"
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9)
  )

save_plot(
  structured_bd_plot,
  file.path(supplement_dir, "GBM_author_label_revised_BD_vs_structured_fraction_supp"),
  width = 8,
  height = 4.2
)


################################################################################
# Supporting QC plot: author label coverage
################################################################################

coverage$has_author_labels <- coverage$n_labeled_spots > 0

coverage_plot <- ggplot(
  coverage,
  aes(x = sample_id, y = prop_labeled_spots, fill = has_author_labels)
) +
  geom_col(width = 0.75) +
  coord_flip() +
  facet_wrap(~ author_label_display, ncol = 1) +
  scale_fill_manual(values = c("FALSE" = "gray75", "TRUE" = "#1B9E77"), guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    x = "Sample",
    y = "Spots with author labels"
  )

save_plot(
  coverage_plot,
  file.path(supplement_dir, "GBM_author_label_revised_coverage_supp"),
  width = 6.5,
  height = 8
)


################################################################################
# Supplementary Figure
################################################################################

summary_row_plot <- wrap_plots(
  relative_bd_label_plot,
  structured_bd_plot,
  nrow = 1,
  guides = "collect"
)

s5_plot <- wrap_plots(
  map_pair_plot,
  summary_row_plot,
  ncol = 1,
  heights = c(1.35, 1),
  guides = "collect"
) +
  plot_layout(heights = c(1.35, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

save_plot(
  s5_plot,
  file.path(supplement_dir, "GBM_author_label_biology_supp"),
  width = 12,
  height = 10
)
