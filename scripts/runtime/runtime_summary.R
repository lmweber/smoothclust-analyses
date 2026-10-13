################################################################################
# Summary of runtimes for GBM and de Oliveira datasets
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(ggplot2)
library(patchwork)

results_dir <- "results"
out_dir <- file.path(results_dir, "runtime")
plots_dir <- file.path("plots", "runtime")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

gbm <- read.csv(
  file.path(results_dir, "GBM", "GBM_runtime_core_smoothclust_summary.csv")
)
external <- read.csv(
  file.path(results_dir, "GBM", "GBM_runtime_external_clustering_summary.csv")
)
hd <- read.csv(
  file.path(results_dir, "deOliveira", "deOliveira_VisiumHD_runtime.csv")
)
xe <- read.csv(
  file.path(results_dir, "deOliveira", "deOliveira_Xenium_runtime.csv")
)

stages <- c("smoothclust", "feature reduction", "clustering and BD")
high_resolution_runtime <- function(x, workflow, platform) {
  data.frame(
    workflow = workflow,
    platform = platform,
    n_sections = 1,
    n_locations = x$n_locations[1],
    n_genes = x$n_genes[1],
    n_threads = x$n_threads[1],
    stage = stages,
    seconds = c(
      x$seconds[x$stage == "smoothclust"],
      x$seconds[x$stage == "sketched_pca_and_full_projection"],
      x$seconds[x$stage == "sketched_kmeans_and_full_assignment"] +
        x$seconds[x$stage == "boundary_density_smoothclust"]
    )
  )
}

# GBM counts shared steps once per section for the selected clustering setting.
# GBM feature reduction includes HVG selection; HD/Xenium genes were preselected.
# Core totals exclude import, output, plotting, author-label BD, and diagnostics.
runtime <- rbind(
  data.frame(
    workflow = sprintf("GBM Visium (%d sections)", gbm$n_samples_after_filter),
    platform = "Visium",
    n_sections = gbm$n_samples_after_filter,
    n_locations = gbm$n_spots_after_filter,
    n_genes = gbm$n_genes_after_filter,
    n_threads = gbm$requested_threads,
    stage = stages,
    seconds = c(
      gbm$comparable_single_k_smoothclust_seconds,
      gbm$comparable_single_k_preprocessing_seconds,
      gbm$comparable_single_k_clustering_bd_seconds
    )
  ),
  high_resolution_runtime(hd, "P2CRC Visium HD", "Visium HD"),
  high_resolution_runtime(xe, "P2CRC Xenium", "Xenium")
)

runtime$total_seconds <- ave(runtime$seconds, runtime$workflow, FUN = sum)
runtime$timing_scope <- ifelse(
  runtime$platform == "Visium",
  sprintf("sum of per-section K=%d workflow components", gbm$comparable_n_clusters),
  paste(
    "full-location smoothing, 30000-location PCA/k-means fit,",
    "full projection/assignment, and full-location BD"
  )
)
write.csv(
  runtime,
  file.path(out_dir, "runtime_core_workflow_components.csv"),
  row.names = FALSE
)

runtime_summary <- unique(
  runtime[, setdiff(names(runtime), c("stage", "seconds"))]
)
write.csv(
  runtime_summary,
  file.path(out_dir, "runtime_core_workflow_summary.csv"),
  row.names = FALSE
)

runtime_metadata <- data.frame(
  field = c(
    "timing_design",
    "source_runs_per_workflow", "peak_memory", "thread_scope",
    "timing_hardware_source", "timing_hardware_model", "timing_processor",
    "timing_processor_cores", "timing_installed_memory"
  ),
  value = c(
    paste0(
      "GBM: component-derived K=", gbm$comparable_n_clusters,
      " total from one K-grid run; HD/Xenium: one observed run each"
    ),
    "1",
    "not measured",
    "Thread counts refer to nearest-neighbor search.",
    "manually specified for the recorded runs",
    "2023 MacBook Pro",
    "Apple M3 Max",
    "14 cores",
    "36 GB"
  )
)
write.csv(
  runtime_metadata,
  file.path(out_dir, "runtime_reproducibility_metadata.csv"),
  row.names = FALSE
)

# Stage levels define the left-to-right order of each stacked bar.
runtime$stage <- factor(runtime$stage, levels = stages)
runtime$workflow_label <- paste0(
  runtime$workflow, "\n",
  format(
    runtime$n_locations, big.mark = ",", scientific = FALSE, trim = TRUE
  ),
  " locations; ",
  format(runtime$n_genes, big.mark = ",", scientific = FALSE, trim = TRUE),
  " genes\n", runtime$n_threads, " NN threads; ",
  ifelse(
    runtime$platform == "Visium",
    paste0("per-section K=", gbm$comparable_n_clusters, " fit"),
    "30,000-location fit"
  )
)
# Order by total runtime; coord_flip places the longest workflow at the top.
label_order <- unique(runtime$workflow_label[order(runtime$total_seconds)])
runtime$workflow_label <- factor(runtime$workflow_label, levels = label_order)

runtime_colors <- c(
  "smoothclust" = "#0072B2",
  "feature reduction" = "#56B4E9",
  "clustering and BD" = "#009E73"
)

p_core <- ggplot(runtime, aes(x = workflow_label, y = seconds, fill = stage)) +
  geom_col(width = 0.62, position = position_stack(reverse = TRUE)) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = runtime_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(
    x = NULL,
    y = "Elapsed time (seconds)",
    fill = NULL,
    title = "Core smoothclust workflow components",
    subtitle = "Descriptive timings; workloads, requested threads, and fitting scope differ"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )

# Reuse recorded comparison totals; no method or timing code is run here.
method_runtime <- rbind(
  data.frame(
    method = "smoothclust",
    runtime_seconds = gbm$comparable_single_k_elapsed_seconds,
    setting = paste0("K = ", gbm$comparable_n_clusters)
  ),
  data.frame(
    method = external$method,
    runtime_seconds = external$total_elapsed_seconds,
    setting = ifelse(
      is.na(external$n_clusters_requested), "default",
      paste0("K = ", external$n_clusters_requested)
    )
  )
)
method_names <- c(
  smoothclust = "smoothclust", bayesspace = "BayesSpace",
  banksy = "BANKSY", clustsignal = "ClustSIGNAL"
)
method_runtime$method_display <- factor(
  method_names[method_runtime$method], levels = unname(method_names)
)
method_runtime$method_label <- factor(
  paste(method_runtime$method_display, method_runtime$setting, sep = "\n"),
  levels = paste(
    method_names, method_runtime$setting[match(names(method_names), method_runtime$method)],
    sep = "\n"
  )
)
method_runtime$runtime_display <- ifelse(
  method_runtime$runtime_seconds < 60,
  sprintf("%.1f s", method_runtime$runtime_seconds),
  ifelse(
    method_runtime$runtime_seconds < 3600,
    sprintf("%.1f min", method_runtime$runtime_seconds / 60),
    sprintf("%.2f h", method_runtime$runtime_seconds / 3600)
  )
)
method_palette <- c(
  "smoothclust" = "#0072B2", "BayesSpace" = "#D55E00",
  "BANKSY" = "#009E73", "ClustSIGNAL" = "#CC79A7"
)
p_method <- ggplot(
  method_runtime,
  aes(x = method_label, y = runtime_seconds, colour = method_display)
) +
  geom_point(size = 3, show.legend = FALSE) +
  geom_text(
    aes(label = runtime_display), colour = "black", vjust = -0.8,
    size = 3.2, show.legend = FALSE
  ) +
  scale_colour_manual(values = method_palette, drop = FALSE) +
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000),
    expand = expansion(mult = c(0.08, 0.18))
  ) +
  labs(
    title = "GBM runtime by clustering method",
    subtitle = paste(gbm$n_samples_after_filter, "section workflows; method timings include BD"),
    x = NULL, y = "Elapsed time (seconds, log10 scale)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

runtime_figure <- (p_core / p_method) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(
    tag_levels = "A",
    caption = paste0(
      "A: GBM component-derived K=", gbm$comparable_n_clusters, " totals over ",
      gbm$n_samples_after_filter, " independent sections; HVG selection and PCA included.\n",
      "HD/Xenium: one observed run each; 30,000-location fit plus full projection/assignment.\n",
      "Core totals exclude import, output, plotting, diagnostics, and HD author-label BD.\n",
      "B: smoothclust reuses the component-derived GBM total; external totals sum per-section timers.\n",
      "Requested K is shown for each method; ClustSIGNAL uses defaults.\n",
      "Method timings include method-specific preprocessing/reduction, clustering/labeling, and BD;\n",
      "exclude shared filtering, import, output, and plotting. GBM requested ",
      paste(unique(external$requested_threads), collapse = ", "),
      " threads where supported.\n",
      "Hardware: 2023 MacBook Pro, Apple M3 Max (14 cores), 36 GB memory."
    ),
    theme = theme(
      plot.tag = element_text(face = "bold"),
      plot.caption = element_text(size = 8, hjust = 0, lineheight = 1.1)
    )
  )

ggsave(
  file.path(plots_dir, "runtime_core_workflow_summary.png"),
  runtime_figure, width = 11, height = 9.5, dpi = 300
)
ggsave(
  file.path(plots_dir, "runtime_core_workflow_summary.pdf"),
  runtime_figure, width = 11, height = 9.5
)
