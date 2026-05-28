################################################################################
# Script to plot runtime summaries
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(ggplot2)
library(patchwork)

dataset_id <- "GBM"

repo_dir <- "."
results_dir <- file.path(repo_dir, "results", dataset_id)
plots_dir <- file.path(repo_dir, "plots", dataset_id)
runtime_plots_dir <- file.path(plots_dir, "runtime")
summary_plots_dir <- file.path(runtime_plots_dir, "summary_plots")
figure_dir <- file.path(runtime_plots_dir, "figure")

dir.create(summary_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

core_runtime_summary_file <- file.path(
  results_dir,
  "GBM_runtime_core_smoothclust_summary.csv"
)
external_runtime_summary_file <- file.path(
  results_dir,
  "GBM_runtime_external_clustering_summary.csv"
)

supplementary_table_file <- file.path(
  results_dir,
  "GBM_runtime_supplementary_table.csv"
)
method_comparison_file <- file.path(
  results_dir,
  "GBM_runtime_method_comparison.csv"
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

format_runtime <- function(seconds) {
  ifelse(
    is.na(seconds),
    NA_character_,
    ifelse(
      seconds < 60,
      sprintf("%.1f s", seconds),
      ifelse(
        seconds < 3600,
        sprintf("%.1f min", seconds / 60),
        sprintf("%.2f h", seconds / 3600)
      )
    )
  )
}

format_integer <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}


# load data

input_files <- c(core_runtime_summary_file, external_runtime_summary_file)

core_runtime <- read.csv(core_runtime_summary_file, stringsAsFactors = FALSE)
external_runtime <- read.csv(external_runtime_summary_file, stringsAsFactors = FALSE)


################################################################################
# Supplementary table data
################################################################################

core_runtime_seconds <- core_runtime$comparable_single_k_elapsed_seconds

runtime_table <- data.frame(
  dataset_id = core_runtime$dataset_id,
  n_samples = core_runtime$n_samples_after_filter,
  n_spots_processed = core_runtime$n_spots_input,
  n_genes_processed = core_runtime$n_genes_input,
  n_spots_analysis = core_runtime$n_spots_after_filter,
  n_genes_analysis = core_runtime$n_genes_after_filter,
  analysis = "smoothclust",
  n_clusters = core_runtime$comparable_n_clusters,
  smoothing_method = core_runtime$smoothing_method,
  smoothing_bandwidth = core_runtime$smoothing_bandwidth,
  boundary_k = core_runtime$boundary_k,
  requested_threads = core_runtime$requested_threads,
  runtime_seconds = core_runtime_seconds,
  runtime_minutes = core_runtime_seconds / 60,
  runtime_display = format_runtime(core_runtime_seconds),
  stringsAsFactors = FALSE
)

write.csv(
  runtime_table,
  file = supplementary_table_file,
  row.names = FALSE
)


################################################################################
# Method-comparison data
################################################################################

method_runtime <- rbind(
  data.frame(
    dataset_id = core_runtime$dataset_id,
    method = "smoothclust",
    method_display = "smoothclust",
    runtime_seconds = core_runtime_seconds,
    n_clusters_requested = paste0("K = ", core_runtime$comparable_n_clusters),
    timing_scope = "single K, derived from core workflow timers",
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset_id = external_runtime$dataset_id,
    method = external_runtime$method,
    method_display = external_runtime$method_display,
    runtime_seconds = external_runtime$total_elapsed_seconds,
    n_clusters_requested = ifelse(
      is.na(external_runtime$n_clusters_requested),
      "default",
      paste0("K = ", external_runtime$n_clusters_requested)
    ),
    timing_scope = "external clustering plus BD calculation",
    stringsAsFactors = FALSE
  )
)

method_runtime$runtime_minutes <- method_runtime$runtime_seconds / 60
method_runtime$runtime_display <- format_runtime(method_runtime$runtime_seconds)

method_levels <- c("smoothclust", "BayesSpace", "Banksy", "ClustSIGNAL")
method_runtime$method_display <- factor(
  method_runtime$method_display,
  levels = method_levels
)
method_runtime <- method_runtime[order(method_runtime$method_display), ]

write.csv(
  method_runtime,
  file = method_comparison_file,
  row.names = FALSE
)


################################################################################
# Plot runtime summaries
################################################################################

theme_set(
  theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
)

runtime_palette <- c(
  "smoothclust" = "#0072B2",
  "BayesSpace" = "#D55E00",
  "Banksy" = "#009E73",
  "ClustSIGNAL" = "#CC79A7"
)

dataset_runtime_plot <- ggplot(
  runtime_table,
  aes(x = reorder(dataset_id, runtime_minutes), y = runtime_minutes)
) +
  geom_col(width = 0.55, fill = runtime_palette[["smoothclust"]]) +
  geom_text(
    aes(label = runtime_display),
    hjust = -0.08,
    size = 3.2
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Core smoothclust runtime",
    subtitle = paste0(
      "Single-K runtime, K = ",
      runtime_table$n_clusters,
      "; ",
      format_integer(runtime_table$n_samples),
      " samples, ",
      format_integer(runtime_table$n_spots_processed),
      " spots"
    ),
    x = NULL,
    y = "Runtime (minutes)"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  )

method_runtime_plot <- ggplot(
  method_runtime,
  aes(x = method_display, y = runtime_minutes, fill = method_display)
) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(
    aes(label = runtime_display),
    vjust = -0.35,
    size = 3.0
  ) +
  scale_fill_manual(values = runtime_palette, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Runtime by clustering method",
    subtitle = "GBM sections; external methods timed with BD calculation",
    x = NULL,
    y = "Runtime (minutes)"
  ) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.major.x = element_blank()
  )

runtime_figure <- dataset_runtime_plot + method_runtime_plot +
  plot_layout(widths = c(0.9, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

save_plot(
  runtime_figure,
  file.path(figure_dir, "GBM_supp_runtime_summary"),
  width = 10,
  height = 4.2
)
