################################################################################
# Script to compare runtimes by clustering method
# Greenwald et al. (2024) GBM/glioma Visium dataset
################################################################################

# run from repository root directory

library(ggplot2)

results_dir <- file.path("results", "GBM")
plots_dir <- file.path("plots", "GBM", "runtime_comparison")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)


# helper function

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


################################################################################
# Method-comparison data
################################################################################

core_runtime <- read.csv(
  file.path(results_dir, "GBM_runtime_core_smoothclust_summary.csv")
)
external_runtime <- read.csv(
  file.path(results_dir, "GBM_runtime_external_clustering_summary.csv")
)

method_runtime <- rbind(
  data.frame(
    dataset_id = core_runtime$dataset_id,
    method = "smoothclust",
    method_display = "smoothclust",
    runtime_seconds = core_runtime$comparable_single_k_elapsed_seconds,
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
method_runtime$method_display <- factor(
  method_runtime$method_display,
  levels = c("smoothclust", "BayesSpace", "Banksy", "ClustSIGNAL")
)
method_runtime <- method_runtime[order(method_runtime$method_display), ]

write.csv(
  method_runtime,
  file.path(results_dir, "GBM_runtime_method_comparison.csv"),
  row.names = FALSE
)


################################################################################
# Plot runtimes
################################################################################

runtime_palette <- c(
  "smoothclust" = "#0072B2",
  "BayesSpace" = "#D55E00",
  "Banksy" = "#009E73",
  "ClustSIGNAL" = "#CC79A7"
)

# Points avoid an arbitrary bar baseline on the logarithmic scale.
method_runtime_plot <- ggplot(
  method_runtime,
  aes(x = method_display, y = runtime_seconds, colour = method_display)
) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text(
    aes(label = runtime_display),
    colour = "black", vjust = -0.8, size = 4, show.legend = FALSE
  ) +
  scale_colour_manual(values = runtime_palette, drop = FALSE) +
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000),
    expand = expansion(mult = c(0.08, 0.18))
  ) +
  labs(
    title = "Runtime by clustering method",
    subtitle = "GBM sections; external methods timed with BD calculation",
    x = NULL,
    y = "Runtime (seconds, log10 scale)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.major.x = element_blank()
  )

ggsave(
  file.path(plots_dir, "GBM_supp_runtime_summary.png"),
  plot = method_runtime_plot, width = 6.5, height = 4.2, dpi = 300
)
ggsave(
  file.path(plots_dir, "GBM_supp_runtime_summary.pdf"),
  plot = method_runtime_plot, width = 6.5, height = 4.2
)
