################################################################################
# Plot Kuppe/Miyoshi boundary-density neighborhood sensitivity
################################################################################

# Run from the smoothclust-analyses repository root, after calculation/checks.

library(ggplot2)
library(patchwork)

results_dir <- file.path("results", "BD_neighborhood_sensitivity")
plot_dir <- file.path("plots", "BD_neighborhood_sensitivity")
bd <- readRDS(file.path(results_dir, "BD_k_sensitivity_section_summary.rds"))
comparisons <- readRDS(file.path(
  results_dir, "BD_k_sensitivity_comparison_summary.rds"
))

group_levels <- c(
  "myogenic-enriched", "ischaemic-enriched", "fibrotic-enriched"
)
group_labels <- c(
  "Myogenic-enriched", "Ischaemic-enriched", "Fibrotic-enriched"
)
group_palette <- setNames(c("#0072B2", "#D55E00", "#009E73"), group_levels)
theme_set(theme_bw(base_size = 14))
panels <- list()
omission_notes <- character()

for (dataset_id in c("Kuppe", "Miyoshi")) {
  dataset_bd <- bd[bd$dataset_id == dataset_id, ]
  baseline <- dataset_bd[dataset_bd$boundary_k == 6, ]
  finite_values <- dataset_bd$relative_boundary_density[
    is.finite(dataset_bd$relative_boundary_density)
  ]
  if (!length(finite_values)) {
    stop(dataset_id, ": no defined relative BD values to plot.")
  }
  row_limits <- range(finite_values)
  padding <- if (diff(row_limits) > 0) 0.05 * diff(row_limits) else 0.05
  row_limits <- row_limits + c(-padding, padding)

  for (k in c(18, 36)) {
    larger <- dataset_bd[dataset_bd$boundary_k == k, ]
    # Calculation saves sections in the same order for each k.
    stopifnot(identical(larger$section_id, baseline$section_id))
    paired <- data.frame(
      section_id = baseline$section_id,
      sample_group = baseline$sample_group,
      relative_bd_k6 = baseline$relative_boundary_density,
      relative_bd_larger_k = larger$relative_boundary_density
    )
    defined <- is.finite(paired$relative_bd_k6) &
      is.finite(paired$relative_bd_larger_k)
    comparison <- comparisons[
      comparisons$dataset_id == dataset_id & comparisons$boundary_k == k &
        comparisons$metric == "relative_boundary_density",
    ]
    if (any(!defined)) {
      message(dataset_id, " k = ", k, ": undefined relative BD; omitted sections: ",
              paste(paired$section_id[!defined], collapse = ", "))
      omission_notes <- c(omission_notes, sprintf(
        "%s k = %d: %d/%d omitted (%s)", dataset_id, k,
        sum(!defined), nrow(paired),
        paste(paired$section_id[!defined], collapse = ", ")
      ))
    }
    rho <- comparison$spearman_rank_correlation
    rho_label <- if (is.finite(rho)) sprintf("%.3f", rho) else "undefined"
    if (!is.finite(rho)) {
      message(dataset_id, " k = ", k, ": Spearman correlation is undefined.")
    }
    p <- ggplot(
      paired[defined, ], aes(relative_bd_k6, relative_bd_larger_k)
    ) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
      coord_equal(xlim = row_limits, ylim = row_limits, expand = FALSE) +
      labs(
        title = dataset_id,
        subtitle = sprintf("Spearman rho = %s; n = %d", rho_label, sum(defined)),
        x = "Relative boundary density (k = 6)",
        y = sprintf("Relative boundary density (k = %d)", k)
      ) +
      theme(plot.title = element_text(face = "bold"))
    if (dataset_id == "Kuppe") {
      p <- p + geom_point(aes(color = sample_group), size = 2.4, alpha = 0.9) +
        scale_color_manual(
          values = group_palette, limits = group_levels, labels = group_labels,
          drop = FALSE
        ) +
        labs(color = "Histomorphological sample category") +
        guides(color = guide_legend(nrow = 1))
    } else {
      p <- p + geom_point(size = 2.4, color = "grey30", alpha = 0.9)
    }
    panels[[length(panels) + 1L]] <- p
  }
}

caption_lines <- character()
if (length(omission_notes)) {
  caption_lines <- strwrap(
    paste("Undefined relative BD:", paste(omission_notes, collapse = "; ")),
    width = 95
  )
}
caption <- if (length(caption_lines)) paste(caption_lines, collapse = "\n") else NULL
figure <- wrap_plots(panels, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Boundary density neighborhood sensitivity",
    caption = caption, tag_levels = "A"
  ) &
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
figure_path <- file.path(plot_dir, "Kuppe_Miyoshi_BD_k_sensitivity")
figure_height <- 9 + 0.16 * max(0, length(caption_lines) - 1)
ggsave(paste0(figure_path, ".pdf"), figure, width = 10, height = figure_height,
       bg = "white")
ggsave(paste0(figure_path, ".png"), figure, width = 10, height = figure_height,
       dpi = 300, bg = "white")
