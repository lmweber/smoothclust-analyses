################################################################################
# Simulation: oversmoothing
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(smoothclust)
library(ggplot2)
library(patchwork)

set.seed(123)

results_dir <- file.path("results", "sim_oversmoothing")
plots_dir <- file.path("plots", "sim_oversmoothing")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

grid_size <- 60L
n_genes <- 120L
n_markers <- 15L
n_hvgs <- 75L
n_pcs <- 10L
n_clusters <- 3L
boundary_k <- 4L
bandwidths <- c(0, 0.025, 0.04, 0.06, 0.085, 0.12, 0.17, 0.23, 0.30)
base_mean_gamma_shape <- 1.5
base_mean_gamma_rate <- 1.5
library_factor_log_sd <- 0.35
negative_binomial_size <- 1
marker_fold_change <- 2

coords <- as.matrix(expand.grid(x = seq_len(grid_size), y = seq_len(grid_size)))
truth <- ifelse(coords[, "x"] <= grid_size / 2, "A", "B")
island <- (coords[, "x"] - 16)^2 + (coords[, "y"] - 30)^2 <= 5^2
truth[island] <- "C"
truth <- factor(truth, levels = c("A", "B", "C"))

base_mean <- rgamma(
  n_genes, shape = base_mean_gamma_shape, rate = base_mean_gamma_rate
)
profile <- matrix(base_mean, nrow = n_genes, ncol = n_clusters)
for (j in seq_len(n_clusters)) {
  marker_ix <- seq.int((j - 1L) * n_markers + 1L, j * n_markers)
  profile[marker_ix, j] <- profile[marker_ix, j] * marker_fold_change
}
library_factors <- exp(rnorm(nrow(coords), sd = library_factor_log_sd))
library_factors <- library_factors / mean(library_factors)
mu <- profile[, as.integer(truth), drop = FALSE]
mu <- sweep(mu, 2, library_factors, "*")
counts <- matrix(
  rnbinom(
    length(mu), mu = as.vector(mu), size = negative_binomial_size
  ),
  nrow = n_genes,
  dimnames = list(paste0("gene", seq_len(n_genes)),
                  paste0("location", seq_len(nrow(coords))))
)
size_factors <- colSums(counts)
size_factors <- size_factors / mean(size_factors)
logcounts <- log2(sweep(counts, 2, size_factors, "/") + 1)
hvg_ix <- order(apply(logcounts, 1, var), decreasing = TRUE)[seq_len(n_hvgs)]

align_clusters <- function(clusters, truth) {
  permutations <- rbind(
    c("A", "B", "C"), c("A", "C", "B"), c("B", "A", "C"),
    c("B", "C", "A"), c("C", "A", "B"), c("C", "B", "A")
  )
  scores <- apply(
    permutations, 1,
    function(x) sum(x[clusters] == as.character(truth))
  )
  factor(permutations[which.max(scores), clusters], levels = c("A", "B", "C"))
}

small_domain_f1 <- function(clusters, truth) {
  target <- truth == "C"
  max(vapply(unique(clusters), function(x) {
    called <- clusters == x
    2 * sum(target & called) / (sum(target) + sum(called))
  }, numeric(1)))
}

truth_bd <- boundary_density(
  coords, truth, k = boundary_k, adjust = TRUE
)
metrics <- vector("list", length(bandwidths))
map_data <- vector("list", length(bandwidths))

for (i in seq_along(bandwidths)) {
  bandwidth <- bandwidths[i]
  expression <- if (bandwidth == 0) {
    logcounts
  } else {
    smoothclust(
      counts,
      spatial_coords = coords,
      method = "uniform",
      bandwidth = bandwidth,
      normalize = TRUE,
      log = TRUE,
      size_factors = size_factors
    )
  }
  pca <- prcomp(
    t(as.matrix(expression[hvg_ix, , drop = FALSE])),
    center = TRUE, scale. = TRUE, rank. = n_pcs
  )$x[, seq_len(n_pcs), drop = FALSE]
  set.seed(123)
  clusters <- kmeans(
    pca, centers = n_clusters, nstart = 50, iter.max = 200,
    algorithm = "MacQueen"
  )$cluster
  inferred <- align_clusters(clusters, truth)
  bd <- boundary_density(
    coords, clusters, k = boundary_k, adjust = TRUE
  )

  metrics[[i]] <- data.frame(
    bandwidth = bandwidth,
    radius_grid_units = bandwidth * (grid_size - 1),
    adjusted_rand_index = mclust::adjustedRandIndex(truth, clusters),
    small_domain_f1 = small_domain_f1(clusters, truth),
    boundary_density = bd$boundary_density,
    expected_boundary_density = bd$expected_boundary_density,
    relative_boundary_density = bd$relative_boundary_density,
    smallest_cluster_n = min(table(clusters))
  )
  map_data[[i]] <- data.frame(
    x = coords[, "x"], y = coords[, "y"], label = inferred,
    bandwidth = bandwidth
  )
}

metrics <- do.call(rbind, metrics)
maps <- do.call(rbind, map_data)
write.csv(
  metrics,
  file.path(results_dir, "simulation_oversmoothing_metrics.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    seed = 123L, grid_size = grid_size, n_locations = nrow(coords),
    n_genes = n_genes, n_markers_per_domain = n_markers,
    n_hvgs = n_hvgs, n_pcs = n_pcs, n_clusters = n_clusters,
    boundary_k = boundary_k, island_n = sum(truth == "C"),
    base_mean_gamma_shape = base_mean_gamma_shape,
    base_mean_gamma_rate = base_mean_gamma_rate,
    library_factor_log_sd = library_factor_log_sd,
    library_factors_mean_centered = TRUE,
    normalization_size_factors = "observed library sums",
    negative_binomial_size = negative_binomial_size,
    marker_fold_change = marker_fold_change
  ),
  file.path(results_dir, "simulation_oversmoothing_parameters.csv"),
  row.names = FALSE
)

domain_colors <- c(A = "#0072B2", B = "#D55E00", C = "#009E73")
theme_map <- theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10)
  )

best_ix <- which.max(metrics$adjusted_rand_index)
representative <- unique(c(0, metrics$bandwidth[best_ix], max(bandwidths)))
truth_df <- data.frame(
  x = coords[, "x"], y = coords[, "y"], label = truth
)
p_truth <- ggplot(truth_df, aes(x, y, fill = label)) +
  geom_raster() +
  scale_fill_manual(values = domain_colors, name = "Domain") +
  guides(fill = guide_legend(ncol = 1, title.position = "top")) +
  coord_fixed() +
  labs(title = "True domains", subtitle = sprintf(
    "Broad A/B domains\n%d-location C island", sum(truth == "C")
  )) +
  theme_map

map_plots <- lapply(representative, function(bandwidth) {
  result <- metrics[metrics$bandwidth == bandwidth, ]
  ggplot(
    maps[maps$bandwidth == bandwidth, ],
    aes(x, y, fill = label)
  ) +
    geom_raster() +
    scale_fill_manual(values = domain_colors, drop = FALSE) +
    coord_fixed() +
    guides(fill = "none") +
    labs(
      title = if (bandwidth == 0) {
        "No smoothing"
      } else {
        sprintf("Bandwidth = %.1f%%", 100 * bandwidth)
      },
      subtitle = sprintf(
        "ARI %.2f\nIsland F1 %.2f", result$adjusted_rand_index,
        result$small_domain_f1
      )
    ) +
    theme_map
})

metric_long <- rbind(
  data.frame(
    bandwidth = metrics$bandwidth, metric = "ARI",
    value = metrics$adjusted_rand_index
  ),
  data.frame(
    bandwidth = metrics$bandwidth, metric = "Island F1",
    value = metrics$small_domain_f1
  )
)
p_recovery <- ggplot(
  metric_long, aes(100 * bandwidth, value, color = metric, shape = metric)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_color_manual(values = c(ARI = "#0072B2", "Island F1" = "#009E73")) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = "Domain recovery",
    x = "Smoothing bandwidth (% of section width)", y = "Recovery",
    color = NULL, shape = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

p_raw <- ggplot(metrics, aes(100 * bandwidth, boundary_density)) +
  geom_hline(
    yintercept = truth_bd$boundary_density, color = "#666666",
    linetype = "dashed"
  ) +
  geom_line(color = "#D55E00", linewidth = 0.7) +
  geom_point(color = "#D55E00", size = 2) +
  scale_y_log10() +
  labs(
    title = "Inferred boundary density",
    subtitle = "Dashed line: true labels",
    x = "Smoothing bandwidth (% of section width)", y = "Raw BD (log scale)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_relative <- ggplot(metrics, aes(100 * bandwidth, relative_boundary_density)) +
  geom_hline(yintercept = 1, color = "#999999", linetype = "dotted") +
  geom_hline(
    yintercept = truth_bd$relative_boundary_density, color = "#666666",
    linetype = "dashed"
  ) +
  geom_line(color = "#CC79A7", linewidth = 0.7) +
  geom_point(color = "#CC79A7", size = 2) +
  scale_y_log10() +
  labs(
    title = "Relative BD:\ncomposition-adjusted summary",
    subtitle = "Dashed: truth\nDotted: random-label expectation",
    x = "Smoothing bandwidth (% of section width)",
    y = "Relative BD (log scale)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

maps_plot <- (p_truth + wrap_plots(map_plots, nrow = 1) +
  plot_layout(guides = "collect")) &
  theme(legend.position = "right", legend.direction = "vertical")

figure <- maps_plot /
  (p_recovery + p_raw + p_relative) +
  plot_annotation(
    title = "Illustrative over-smoothing: lower BD with poorer recovery",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 15),
      plot.tag = element_text(size = 15)
    )
  )

ggsave(
  file.path(plots_dir, "sim_oversmoothing.png"),
  figure, width = 12, height = 7.6, dpi = 300, bg = "white"
)
ggsave(
  file.path(plots_dir, "sim_oversmoothing.pdf"),
  figure, width = 12, height = 7.6, bg = "white"
)
