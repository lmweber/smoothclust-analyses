################################################################################
# Simulation: illustrate composition adjustment underlying relative BD
################################################################################

# Run from the smoothclust-analyses repository root directory.
# Recreates the simulation results and figures in one complete run.

library(smoothclust)
library(ggplot2)
library(patchwork)

results_dir <- file.path("results", "sim_composition")
plots_dir <- file.path("plots", "sim_composition")
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
    adjust = TRUE
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

# Random-label reference: randomize labels while retaining exact counts.
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
      labels, "random_label_reference", sprintf("%.1f", proportion), j,
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

# Composition series: move one straight interface without changing its length.
interface_x <- round(grid_side * null_proportions)
interface_graph_seed <- seed + 6001L
interface_labels <- lapply(interface_x, function(x) {
  ifelse(coords[, "x"] <= x, "A", "B")
})
names(interface_labels) <- paste0(
  round(100 * null_proportions), ":", round(100 * (1 - null_proportions))
)
interface_results <- do.call(rbind, lapply(seq_along(interface_labels), function(i) {
  # Hold the neighbor-search seed fixed across these deterministic label maps.
  set.seed(interface_graph_seed)
  evaluate_labels(
    interface_labels[[i]], "straight_interface_series",
    names(interface_labels)[i], 1L, NA_integer_
  )
}))
interface_results$state <- factor(
  interface_results$state, levels = names(interface_labels)
)
interface_results$interface_x <- interface_x
interface_results$graph_seed <- interface_graph_seed

write.csv(
  null_results,
  file.path(results_dir, "BD_random_label_reference.csv"),
  row.names = FALSE
)
write.csv(
  interface_results,
  file.path(results_dir, "BD_straight_interface_series.csv"),
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
    seed = seed, n_threads = n_threads,
    interface_graph_seed = interface_graph_seed
  ),
  file.path(results_dir, "BD_simulation_parameters.csv"),
  row.names = FALSE
)

# Figure: random-label reference, organization, and composition adjustment.
label_palette <- c(A = "#0072B2", B = "#E69F00")
theme_set(
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      strip.text = element_text(size = 10)
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
      title = as.character(result$state),
      subtitle = sprintf(
        "Raw BD = %.3f\nRelative BD = %.3f",
        result$boundary_density, result$relative_boundary_density
      )
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0.5),
      plot.margin = margin(5.5, 10, 5.5, 10)
    )
})
maps_plot <- wrap_plots(map_plots, nrow = 1)

null_results$composition <- factor(
  paste0(round(100 * null_results$proportion_A), ":", 
         round(100 * (1 - null_results$proportion_A))),
  levels = c("10:90", "30:70", "50:50")
)
null_expectations <- unique(null_results[c(
  "composition", "expected_boundary_density"
)])
null_raw_plot <- ggplot(null_results, aes(composition, boundary_density)) +
  geom_boxplot(
    width = 0.55, outlier.shape = NA, fill = "#D55E00",
    alpha = 0.55, linewidth = 0.35
  ) +
  geom_jitter(width = 0.12, height = 0, size = 0.7, alpha = 0.45) +
  geom_errorbar(
    data = null_expectations,
    aes(x = composition, ymin = expected_boundary_density,
        ymax = expected_boundary_density),
    inherit.aes = FALSE, width = 0.7, linetype = "dashed",
    color = "grey35", linewidth = 0.4
  ) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2,
               fill = "white") +
  labs(
    title = "Random arrangements: raw BD",
    subtitle = "Dashed bars: analytic composition-specific expectation",
    x = "A:B label composition", y = "Raw BD"
  )
null_relative_plot <- ggplot(
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
    title = "Random arrangements: relative BD",
    subtitle = "After composition adjustment; dashed line: reference of 1",
    x = "A:B label composition", y = "Relative BD"
  )
null_plot <- wrap_plots(
  plot_spacer(), null_raw_plot + null_relative_plot, plot_spacer(),
  nrow = 1, widths = c(0.08, 1, 0.08)
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
  scale_color_manual(
    values = c(Observed = "#D55E00", Expected = "grey30"),
    labels = c(Observed = "Raw BD", Expected = "Expected raw BD")
  ) +
  labs(title = "Raw and expected raw BD", x = NULL,
       y = "Raw and expected raw BD", color = NULL) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )
# At fixed label counts, relative BD is a constant scaling of raw BD.
organization_plot <- wrap_plots(
  plot_spacer(), p_raw, plot_spacer(),
  nrow = 1, widths = c(0.08, 1, 0.08)
)

interface_plots <- lapply(seq_along(interface_labels), function(i) {
  result <- interface_results[i, ]
  map <- data.frame(x = coords[, "x"], y = coords[, "y"],
                    label = interface_labels[[i]])
  ggplot(map, aes(x, y, fill = label)) +
    geom_tile() +
    scale_fill_manual(values = label_palette, guide = "none") +
    coord_equal(expand = FALSE) +
    labs(
      title = paste0("A:B = ", result$state),
      subtitle = sprintf(
        "Raw BD = %.3f\nExpected raw BD = %.3f\nRelative BD = %.3f",
        result$boundary_density, result$expected_boundary_density,
        result$relative_boundary_density
      )
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, hjust = 0),
      plot.margin = margin(5.5, 10, 5.5, 10)
    )
})
interface_plot <- wrap_plots(interface_plots, nrow = 1)

tag_group <- function(plot, tag) {
  wrap_elements(full = plot) + labs(tag = tag) +
    theme(
      plot.tag = element_text(face = "bold", size = 14, hjust = 0),
      plot.tag.position = c(0, 1),
      plot.margin = margin(5.5, 5.5, 5.5, 15)
    )
}
figure <- tag_group(maps_plot, "A") / tag_group(null_plot, "B") /
  tag_group(organization_plot, "C") / tag_group(interface_plot, "D") +
  plot_layout(heights = c(0.95, 1.1, 1.25, 1.2))
ggsave(
  file.path(plots_dir, "sim_composition.pdf"),
  figure, width = 11, height = 10.5, bg = "white"
)
ggsave(
  file.path(plots_dir, "sim_composition.png"),
  figure, width = 11, height = 10.5, dpi = 300, bg = "white"
)
