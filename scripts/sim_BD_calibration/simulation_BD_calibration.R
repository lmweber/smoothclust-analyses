################################################################################
# Simulate calibration and spatial-organization behavior of analytic BD
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(smoothclust)
library(ggplot2)
library(patchwork)

results_dir <- file.path("results", "sim_BD_calibration")
plots_dir <- file.path("plots", "sim_BD_calibration")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

grid_side <- 30L
boundary_k <- 4L
n_null_replicates <- 75L
seed <- 123L
n_threads <- 4L

coords <- as.matrix(expand.grid(x = seq_len(grid_side), y = seq_len(grid_side)))
n_spots <- nrow(coords)

evaluate_labels <- function(labels, simulation, state, replicate, label_seed) {
  bd <- boundary_density(
    coords, labels, k = boundary_k, n_threads = n_threads,
    adjust = "analytic"
  )
  counts <- table(factor(labels, levels = c("A", "B")))
  data.frame(
    simulation = simulation, state = state, replicate = replicate,
    proportion_A = counts[["A"]] / n_spots,
    n_A = counts[["A"]], n_B = counts[["B"]], n_spots = n_spots,
    boundary_k = boundary_k, label_seed = label_seed,
    boundary_density = bd$boundary_density,
    expected_boundary_density = bd$expected_boundary_density,
    relative_boundary_density = bd$relative_boundary_density,
    stringsAsFactors = FALSE
  )
}

# Null calibration: independently randomize labels while retaining exact counts.
null_proportions <- c(0.1, 0.3, 0.5)
null_results <- do.call(rbind, lapply(seq_along(null_proportions), function(i) {
  proportion <- null_proportions[i]
  do.call(rbind, lapply(seq_len(n_null_replicates), function(j) {
    label_seed <- seed + i * 10000L + j
    set.seed(label_seed)
    labels <- sample(c(
      rep("A", round(n_spots * proportion)),
      rep("B", n_spots - round(n_spots * proportion))
    ))
    evaluate_labels(
      labels, "null_calibration", sprintf("%.1f", proportion), j,
      label_seed
    )
  }))
}))

# Organization series: change spatial arrangement while retaining 450:450 counts.
coherent <- ifelse(coords[, "x"] <= grid_side / 2, "A", "B")
swap_labels <- function(labels, fraction, swap_seed) {
  set.seed(swap_seed)
  n_each <- round(length(labels) * fraction / 2)
  a <- sample(which(labels == "A"), n_each)
  b <- sample(which(labels == "B"), n_each)
  labels[c(a, b)] <- labels[c(b, a)]
  labels
}

organization_label_seeds <- c(
  NA_integer_, seed + 5001L, seed + 5002L, seed + 5003L, NA_integer_
)
set.seed(organization_label_seeds[4])
random_labels <- sample(coherent)
organization_labels <- list(
  "Coherent" = coherent,
  "10% swapped" = swap_labels(coherent, 0.1, organization_label_seeds[2]),
  "30% swapped" = swap_labels(coherent, 0.3, organization_label_seeds[3]),
  "Random" = random_labels,
  "Checkerboard" = ifelse(rowSums(coords) %% 2 == 0, "A", "B")
)
organization_results <- do.call(rbind, lapply(
  seq_along(organization_labels), function(i) {
    evaluate_labels(
      organization_labels[[i]], "organization_series",
      names(organization_labels)[i], 1L, organization_label_seeds[i]
    )
  }
))
organization_results$state <- factor(
  organization_results$state, levels = names(organization_labels)
)

write.csv(
  null_results,
  file.path(results_dir, "BD_null_calibration.csv"),
  row.names = FALSE
)
write.csv(
  organization_results,
  file.path(results_dir, "BD_organization_series.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    grid_side = grid_side, n_spots = n_spots, boundary_k = boundary_k,
    n_null_replicates = n_null_replicates,
    seed = seed, n_threads = n_threads
  ),
  file.path(results_dir, "BD_simulation_parameters.csv"),
  row.names = FALSE
)

# Figure: spatial configurations, null calibration, and organization response.
label_palette <- c(A = "#0072B2", B = "#E69F00")
theme_set(
  theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      strip.text = element_text(size = 8)
    )
)

map_data <- do.call(rbind, lapply(seq_along(organization_labels), function(i) {
  data.frame(
    x = coords[, "x"], y = coords[, "y"],
    label = organization_labels[[i]],
    state = factor(names(organization_labels)[i],
                   levels = names(organization_labels))
  )
}))
map_plots <- lapply(seq_along(organization_labels), function(i) {
  result <- organization_results[i, ]
  ggplot(map_data[map_data$state == result$state, ], aes(x, y, fill = label)) +
    geom_tile() +
    scale_fill_manual(values = label_palette, guide = "none") +
    coord_equal(expand = FALSE) +
    labs(
      title = paste0(if (i == 1L) "A  " else "", as.character(result$state)),
      subtitle = sprintf(
        "raw BD = %.3f; rBD = %.2f",
        result$boundary_density, result$relative_boundary_density
      )
    ) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
      plot.subtitle = element_text(size = 7.5, hjust = 0.5)
    )
})
maps_plot <- wrap_plots(map_plots, nrow = 1)

null_results$composition <- factor(
  paste0(round(100 * null_results$proportion_A), ":", 
         round(100 * (1 - null_results$proportion_A))),
  levels = c("10:90", "30:70", "50:50")
)
null_plot <- ggplot(
  null_results, aes(composition, relative_boundary_density)
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed", color = "grey35", linewidth = 0.4
  ) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, fill = "#56B4E9",
    alpha = 0.55, linewidth = 0.35
  ) +
  geom_jitter(width = 0.12, height = 0, size = 0.7, alpha = 0.45) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2,
               fill = "white") +
  labs(
    title = "B  Null centering with exact label counts",
    x = "A:B label composition", y = "Relative BD"
  )

p_raw <- ggplot(organization_results, aes(state, group = 1)) +
  geom_line(aes(y = boundary_density, color = "Observed"), linewidth = 0.6) +
  geom_point(aes(y = boundary_density, color = "Observed"), size = 1.8) +
  geom_line(
    aes(y = expected_boundary_density, color = "Expected"),
    linetype = "dashed", linewidth = 0.6
  ) +
  geom_point(
    aes(y = expected_boundary_density, color = "Expected"),
    shape = 1, size = 1.8
  ) +
  scale_color_manual(values = c(Observed = "#D55E00", Expected = "grey30")) +
  labs(title = "C  Raw and expected BD", x = NULL, y = "BD", color = NULL) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )
p_relative <- ggplot(
  organization_results, aes(state, relative_boundary_density, group = 1)
) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey35", linewidth = 0.4) +
  geom_line(color = "#009E73", linewidth = 0.6) +
  geom_point(color = "#009E73", size = 1.8) +
  labs(x = NULL, y = "Relative BD") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
organization_plot <- p_raw + p_relative

figure <- maps_plot / null_plot / organization_plot +
  plot_layout(heights = c(0.95, 1.1, 1.25))
ggsave(
  file.path(plots_dir, "sim_BD_calibration.pdf"),
  figure, width = 11, height = 8.5, bg = "white"
)
ggsave(
  file.path(plots_dir, "sim_BD_calibration.png"),
  figure, width = 11, height = 8.5, dpi = 300, bg = "white"
)
