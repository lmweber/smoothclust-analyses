################################################################################
# Script to plot within-tumor author-label boundary density results
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(ggplot2)
library(patchwork)
library(scales)
library(SpatialExperiment)

dataset_id <- "GBM"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)
plots_dir <- file.path(repo_dir, "plots", dataset_id)
within_tumor_dir <- file.path(
  plots_dir,
  "author_labels",
  "within_tumor"
)

dir.create(within_tumor_dir, recursive = TRUE, showWarnings = FALSE)

section_summary_file <- file.path(
  results_dir,
  "GBM_BD_author_labels_within_tumor_section_summary.csv"
)
composition_file <- file.path(
  results_dir,
  "GBM_BD_author_labels_within_tumor_composition.csv"
)
spe_plot_file <- file.path(
  results_dir,
  "spe_GBM_author_labels_BD_plot.rds"
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

require_columns <- function(x, columns, object_name) {
  missing_columns <- setdiff(columns, names(x))
  if (length(missing_columns) > 0L) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

ordered_levels <- function(observed, preferred) {
  observed <- unique(as.character(observed))
  observed <- observed[!is.na(observed) & nzchar(observed)]
  c(intersect(preferred, observed), sort(setdiff(observed, preferred)))
}


# load data

input_files <- c(section_summary_file, composition_file, spe_plot_file)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing required evaluation outputs: ",
    paste(missing_files, collapse = ", ")
  )
}

section_summary <- read.csv(section_summary_file, stringsAsFactors = FALSE)
composition <- read.csv(composition_file, stringsAsFactors = FALSE)
spe <- readRDS(spe_plot_file)

require_columns(
  section_summary,
  c(
    "sample_id",
    "tumor_group",
    "region_group",
    "section_display",
    "tumor_section_order",
    "fraction_structured",
    "fraction_mes_hyp",
    "mp_relative_boundary_density",
    "layer_relative_boundary_density"
  ),
  "section_summary"
)
require_columns(
  composition,
  c(
    "sample_id",
    "tumor_group",
    "section_display",
    "tumor_section_order",
    "metaprogram",
    "fraction_metaprogram"
  ),
  "composition"
)
require_columns(
  as.data.frame(colData(spe)),
  c(
    "sample_id",
    "author_mp",
    "local_bd_author_mp"
  ),
  "colData(spe)"
)

if (nrow(section_summary) == 0L ||
    anyNA(section_summary$sample_id) ||
    anyDuplicated(section_summary$sample_id)) {
  stop("section_summary must contain one uniquely identified row per section.")
}
section_numeric <- c(
  "tumor_section_order",
  "fraction_structured",
  "fraction_mes_hyp",
  "mp_relative_boundary_density",
  "layer_relative_boundary_density"
)
if (any(!is.finite(as.matrix(section_summary[section_numeric])))) {
  stop("Non-finite plotting values found in section_summary.")
}
composition_key <- paste(
  composition$sample_id,
  composition$metaprogram,
  sep = "___"
)
if (nrow(composition) == 0L ||
    anyNA(composition$sample_id) ||
    anyNA(composition$metaprogram) ||
    anyDuplicated(composition_key) ||
    any(!is.finite(composition$fraction_metaprogram))) {
  stop("composition must contain one finite row per section and metaprogram.")
}

tumor_levels <- ordered_levels(
  section_summary$tumor_group,
  c("ZH1007", "ZH1019", "ZH881", "ZH916")
)
section_summary$tumor_group <- factor(
  section_summary$tumor_group,
  levels = tumor_levels
)
composition$tumor_group <- factor(composition$tumor_group, levels = tumor_levels)

preferred_mp_levels <- c(
  "AC",
  "OPC",
  "NPC",
  "MES",
  "MES.Ast",
  "MES.Hyp",
  "Prolif.Metab",
  "Chromatin.Reg",
  "Vasc",
  "Mac",
  "Inflammatory.Mac",
  "Reactive.Ast",
  "Oligo",
  "Neuron"
)
mp_levels <- ordered_levels(composition$metaprogram, preferred_mp_levels)
if (length(mp_levels) == 0L) {
  stop("No dominant metaprogram labels found in the composition table.")
}
composition$metaprogram <- factor(composition$metaprogram, levels = mp_levels)


# plot settings

theme_set(
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "gray90", color = "gray45")
    )
)

mp_palette <- setNames(
  grDevices::hcl.colors(length(mp_levels), palette = "Dark 3"),
  mp_levels
)

label_source_levels <- c("Dominant metaprogram", "State layer")
label_source_palette <- c(
  "Dominant metaprogram" = "#0072B2",
  "State layer" = "#D55E00"
)
label_source_shapes <- c(
  "Dominant metaprogram" = 21,
  "State layer" = 24
)


################################################################################
# Panel A: dominant metaprogram and local-BD maps for selected ZH881 sections
################################################################################

# These three sections were selected before analysis because Greenwald et al.
# used them to illustrate within-tumor differences in ZH881.
map_sample_ids <- c("zh881t1", "zh8812", "zh8811a")
map_section_summary <- section_summary[
  match(map_sample_ids, section_summary$sample_id),
  ,
  drop = FALSE
]
if (any(is.na(map_section_summary$sample_id))) {
  stop("Could not match all preselected ZH881 map sections.")
}

map_section_summary$map_facet <- paste0(
  map_section_summary$section_display,
  "\nMP relative BD = ",
  sprintf("%.2f", map_section_summary$mp_relative_boundary_density)
)
map_facet_levels <- map_section_summary$map_facet

map_df <- do.call(
  rbind,
  lapply(seq_along(map_sample_ids), function(i) {
    sample_id <- map_sample_ids[i]
    sample_ix <- as.character(colData(spe)$sample_id) == sample_id
    if (!any(sample_ix)) {
      stop("No spots found in spe for map section: ", sample_id)
    }
    map_coords <- spatialCoords(spe)[sample_ix, , drop = FALSE]
    data.frame(
      sample_id = sample_id,
      map_facet = map_section_summary$map_facet[i],
      x = map_coords[, 1],
      y = -map_coords[, 2],
      metaprogram = as.character(colData(spe)$author_mp[sample_ix]),
      local_boundary_density = as.numeric(
        colData(spe)$local_bd_author_mp[sample_ix]
      )
    )
  })
)

map_df$map_facet <- factor(map_df$map_facet, levels = map_facet_levels)
map_df$metaprogram <- factor(map_df$metaprogram, levels = mp_levels)
map_theme <- theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.background = element_rect(fill = "gray90", color = "gray45"),
    strip.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10)
  )
map_base <- ggplot(map_df, aes(x = x, y = y)) +
  geom_point(color = "gray90", size = 0.4) +
  coord_fixed() +
  facet_wrap(~ map_facet, nrow = 1) +
  map_theme

mp_map_plot <- map_base +
  geom_point(
    data = map_df[!is.na(map_df$metaprogram), , drop = FALSE],
    aes(color = metaprogram),
    shape = 16,
    size = 0.52,
    show.legend = TRUE
  ) +
  scale_color_manual(values = mp_palette, limits = mp_levels, drop = FALSE) +
  labs(
    title = "Dominant metaprogram maps",
    x = NULL,
    y = NULL,
    color = "Dominant metaprogram"
  ) +
  theme(
    legend.text = element_text(size = 8.5),
    legend.key.height = grid::unit(0.7, "lines")
  ) +
  guides(
    color = guide_legend(
      ncol = 2,
      byrow = TRUE,
      override.aes = list(size = 2.5)
    )
  )

local_bd_map_plot <- map_base +
  geom_point(
    data = map_df[
      is.finite(map_df$local_boundary_density),
      ,
      drop = FALSE
    ],
    aes(color = local_boundary_density),
    size = 0.52
  ) +
  scale_color_gradientn(
    colors = rev(grDevices::hcl.colors(60, palette = "YlOrRd")),
    limits = c(0, 1),
    oob = scales::squish
  ) +
  labs(
    title = "Local boundary density from dominant metaprogram labels",
    x = NULL,
    y = NULL,
    color = "Local BD"
  ) +
  theme(legend.text = element_text(size = 9))

map_pair_plot <- mp_map_plot / local_bd_map_plot +
  plot_layout(heights = c(1, 1))


################################################################################
# Panel B: dominant-metaprogram composition of all five ZH881 sections
################################################################################

zh881_summary <- section_summary[
  section_summary$tumor_group == "ZH881",
  ,
  drop = FALSE
]
zh881_summary <- zh881_summary[
  order(zh881_summary$tumor_section_order),
  ,
  drop = FALSE
]
zh881_section_levels <- zh881_summary$section_display

zh881_composition <- composition[
  composition$tumor_group == "ZH881",
  ,
  drop = FALSE
]
zh881_composition$section_display <- factor(
  zh881_composition$section_display,
  levels = zh881_section_levels
)

composition_plot <- ggplot(
  zh881_composition,
  aes(
    x = section_display,
    y = fraction_metaprogram,
    fill = metaprogram
  )
) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = mp_palette, drop = FALSE, guide = "none") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "ZH881 metaprogram composition",
    x = "Section",
    y = "Fraction of MP-labeled spots"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )


################################################################################
# ZH881 relative BD and context; retain the context percentages in panel B
################################################################################

zh881_bd_long <- data.frame(
  section_display = rep(zh881_summary$section_display, times = 2),
  label_source = rep(label_source_levels, each = nrow(zh881_summary)),
  relative_boundary_density = c(
    zh881_summary$mp_relative_boundary_density,
    zh881_summary$layer_relative_boundary_density
  )
)
zh881_bd_long$section_display <- factor(
  zh881_bd_long$section_display,
  levels = zh881_section_levels
)
zh881_bd_long$label_source <- factor(
  zh881_bd_long$label_source,
  levels = label_source_levels
)

zh881_bd_plot <- ggplot(
  zh881_bd_long,
  aes(
    x = section_display,
    y = relative_boundary_density,
    color = label_source,
    shape = label_source
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray45") +
  geom_point(
    position = position_dodge(width = 0.42),
    size = 2.5,
    stroke = 0.7
  ) +
  scale_color_manual(values = label_source_palette, drop = FALSE) +
  scale_shape_manual(values = c(16, 17), drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.03))) +
  labs(
    title = "ZH881 spatial organization",
    x = NULL,
    y = "Relative boundary density",
    color = "Author label",
    shape = "Author label"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.margin = margin(t = -3, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE)
  )

zh881_context <- data.frame(
  section_display = rep(zh881_summary$section_display, times = 2),
  context = rep(
    c("Structured zone", "MES-Hyp"),
    each = nrow(zh881_summary)
  ),
  fraction = c(
    zh881_summary$fraction_structured,
    zh881_summary$fraction_mes_hyp
  )
)
zh881_context$section_display <- factor(
  zh881_context$section_display,
  levels = zh881_section_levels
)
zh881_context$context <- factor(
  zh881_context$context,
  levels = c("MES-Hyp", "Structured zone")
)
zh881_context$percent_label <- percent(
  zh881_context$fraction,
  accuracy = 1
)
zh881_context$text_color <- ifelse(
  is.finite(zh881_context$fraction) & zh881_context$fraction >= 0.55,
  "white",
  "gray15"
)

zh881_context_plot <- ggplot(
  zh881_context,
  aes(x = section_display, y = context, fill = fraction)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(
    aes(label = percent_label, color = text_color),
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#2171B5",
    limits = c(0, 1),
    oob = scales::squish,
    guide = "none"
  ) +
  scale_color_identity() +
  labs(x = "Section", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.ticks = element_blank(),
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  )

zh881_composition_context_plot <- (
  composition_plot +
    labs(x = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
) / zh881_context_plot +
  plot_layout(heights = c(2.4, 1))


################################################################################
# Panel C: descriptive within-tumor comparison across all 12 sections
################################################################################

section_summary <- section_summary[
  order(section_summary$tumor_group, section_summary$tumor_section_order),
  ,
  drop = FALSE
]

# Reverse each tumor's levels so its first section is displayed at the top.
section_levels <- unlist(
  lapply(
    tumor_levels,
    function(tumor_group) {
      reverse_ix <- which(as.character(section_summary$tumor_group) == tumor_group)
      rev(section_summary$sample_id[reverse_ix])
    }
  ),
  use.names = FALSE
)
section_summary$sample_id <- factor(
  section_summary$sample_id,
  levels = section_levels
)
section_label_lookup <- setNames(
  section_summary$section_display,
  as.character(section_summary$sample_id)
)

all_bd_long <- section_summary[
  rep(seq_len(nrow(section_summary)), times = 2),
  c("sample_id", "tumor_group", "fraction_structured"),
  drop = FALSE
]
all_bd_long$label_source <- rep(
  label_source_levels,
  each = nrow(section_summary)
)
all_bd_long$relative_boundary_density <- c(
  section_summary$mp_relative_boundary_density,
  section_summary$layer_relative_boundary_density
)
all_bd_long$tumor_group <- factor(all_bd_long$tumor_group, levels = tumor_levels)
all_bd_long$label_source <- factor(
  all_bd_long$label_source,
  levels = label_source_levels
)

all_tumor_plot <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray45") +
  geom_segment(
    data = section_summary,
    aes(
      x = mp_relative_boundary_density,
      xend = layer_relative_boundary_density,
      y = sample_id,
      yend = sample_id
    ),
    color = "gray55",
    linewidth = 0.6
  ) +
  geom_point(
    data = all_bd_long,
    aes(
      x = relative_boundary_density,
      y = sample_id,
      shape = label_source,
      fill = fraction_structured
    ),
    color = "gray15",
    size = 2.8,
    stroke = 0.6
  ) +
  facet_wrap(~ tumor_group, nrow = 1, scales = "free_y") +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_discrete(labels = section_label_lookup) +
  scale_shape_manual(values = label_source_shapes, drop = FALSE) +
  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#2171B5",
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = percent_format(accuracy = 1),
    oob = scales::squish
  ) +
  labs(
    title = "Within-tumor section differences",
    x = "Relative boundary density",
    y = NULL,
    shape = "Author label",
    fill = "Structured-zone fraction"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = -2, unit = "pt"),
    strip.text = element_text(face = "bold")
  ) +
  guides(
    shape = guide_legend(nrow = 1, byrow = TRUE),
    fill = guide_colorbar(
      title.position = "top",
      barwidth = grid::unit(2.8, "cm"),
      barheight = grid::unit(0.25, "cm")
    )
  )


################################################################################
# Supplementary figure
################################################################################

# Keep tags on full-width rows while preserving the B/C content insets.
within_tumor_supp_plot <- wrap_plots(
  A = wrap_elements(full = map_pair_plot),
  B = wrap_elements(full = wrap_plots(
    plot_spacer(), zh881_composition_context_plot, plot_spacer(),
    nrow = 1, widths = c(1, 10, 1)
  )),
  C = wrap_elements(full = wrap_plots(
    plot_spacer(), all_tumor_plot, plot_spacer(),
    nrow = 1, widths = c(1, 10, 1)
  )),
  ncol = 1,
  heights = c(2, 1, 0.9)
) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1),
    plot.tag.position = c(0.05, 0.995)
  )

save_plot(
  within_tumor_supp_plot,
  file.path(
    within_tumor_dir,
    "GBM_author_label_within_tumor_supp"
  ),
  width = 15,
  height = 13
)
