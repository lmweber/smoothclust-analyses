################################################################################
# Plot illustrative de Oliveira Visium HD and Xenium analyses
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(ggplot2)
library(ggrastr)
library(patchwork)

dataset_id <- "deOliveira"
results_dir <- file.path("results", dataset_id)
plots_dir <- file.path("plots", dataset_id)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, filename, width, height) {
  ggsave(paste0(filename, ".png"), plot, width = width, height = height,
         dpi = 300, bg = "white")
  ggsave(paste0(filename, ".pdf"), plot, width = width, height = height,
         bg = "white")
}

smoothclust_colors <- setNames(
  hcl.colors(15, palette = "Dynamic"),
  as.character(seq_len(15))
)
author_colors <- setNames(
  hcl.colors(15, palette = "Dark 3"),
  as.character(seq_len(15))
)

hd <- readRDS(
  file.path(results_dir, "deOliveira_VisiumHD_spot_results.rds")
)
xe <- readRDS(
  file.path(results_dir, "deOliveira_Xenium_spot_results.rds")
)

theme_spatial <- theme_void(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 9)
  )

# Rasterize spatial layers to keep PDF files manageable.
hd_author <- ggplot(hd, aes(x = x, y = y, fill = author_cluster)) +
  rasterise(geom_tile(width = 1, height = 1), dpi = 300) +
  scale_fill_manual(values = author_colors, name = "Cluster") +
  scale_y_reverse() +
  coord_fixed() +
  guides(fill = "none") +
  labs(
    title = "Distributed cluster labels",
    subtitle = "Full P2CRC 8-µm Visium HD section"
  ) +
  theme_spatial

hd_smooth <- ggplot(hd, aes(x = x, y = y, fill = smoothclust_label)) +
  rasterise(geom_tile(width = 1, height = 1), dpi = 300) +
  scale_fill_manual(values = smoothclust_colors, name = "Domain") +
  scale_y_reverse() +
  coord_fixed() +
  guides(fill = "none") +
  labs(
    title = "smoothclust-smoothed domains",
    subtitle = "30,000-bin sketched PCA/k-means fit; 15 domains"
  ) +
  theme_spatial

hd_local <- ggplot(hd, aes(x = x, y = y, fill = local_bd_smoothclust)) +
  rasterise(geom_tile(width = 1, height = 1), dpi = 300) +
  scale_fill_viridis_c(limits = c(0, 1), name = "Local BD") +
  scale_y_reverse() +
  coord_fixed() +
  labs(title = "Local boundary density", subtitle = "Eight-neighbor BD graph") +
  theme_spatial

hd_figure <- hd_author + hd_smooth + hd_local +
  plot_annotation(
    title = "de Oliveira P2CRC Visium HD",
    subtitle = paste(
      "Author-cluster and smoothclust colors are independent;",
      "cluster numbers are not matched"
    ),
    theme = theme(text = element_text(size = 14))
  )

xe_smooth <- ggplot(xe, aes(x = x, y = y, color = smoothclust_label)) +
  geom_point_rast(size = 0.08, raster.dpi = 300) +
  scale_color_manual(values = smoothclust_colors, name = "Domain") +
  scale_y_reverse() +
  coord_fixed() +
  guides(color = "none") +
  labs(
    title = "smoothclust-smoothed domains",
    subtitle = "30,000-cell sketched PCA/k-means fit; 15 domains"
  ) +
  theme_spatial

xe_local <- ggplot(xe, aes(x = x, y = y, color = local_bd_smoothclust)) +
  geom_point_rast(size = 0.08, raster.dpi = 300) +
  scale_color_viridis_c(limits = c(0, 1), name = "Local BD") +
  scale_y_reverse() +
  coord_fixed() +
  labs(title = "Local boundary density", subtitle = "Six-neighbor cell graph") +
  theme_spatial

xe_figure <- xe_smooth + xe_local +
  plot_annotation(
    title = "de Oliveira P2CRC Xenium",
    theme = theme(text = element_text(size = 14))
  )

save_plot(
  hd_figure,
  file.path(plots_dir, "deOliveira_VisiumHD_spatial_summary"),
  15, 5.5
)
save_plot(
  xe_figure,
  file.path(plots_dir, "deOliveira_Xenium_spatial_summary"),
  11, 5.5
)


# Compact descriptive metric and runtime summaries

bd <- rbind(
  read.csv(file.path(results_dir, "deOliveira_VisiumHD_BD_summary.csv")),
  read.csv(file.path(results_dir, "deOliveira_Xenium_BD_summary.csv"))
)

# Numeric panel for the retained spatial figure; keep the historical bar export.
bd_numeric <- bd[order(bd$platform, bd$label_source), ]
label_sources <- c(
  "OSTA distributed clusters" = "distributed labels",
  "smoothclust + sketched PCA/k-means" = "smoothclust labels"
)
stopifnot(nrow(bd_numeric) == 3L, !anyNA(label_sources[bd_numeric$label_source]))
table_values <- list(
  c("Platform / label source", paste(
    bd_numeric$platform, label_sources[bd_numeric$label_source], sep = " / "
  )),
  c("BD neighbors\n(k)", bd_numeric$boundary_k),
  c("Raw BD", sprintf("%.4f", bd_numeric$boundary_density)),
  c("Expected\nraw BD", sprintf("%.4f", bd_numeric$expected_boundary_density)),
  c("Relative BD", sprintf("%.4f", bd_numeric$relative_boundary_density)),
  c("Locations", format(
    bd_numeric$n_locations, big.mark = ",", scientific = FALSE, trim = TRUE
  ))
)
table_cells <- data.frame(
  x = rep(c(0.15, 5.3, 6.65, 8.15, 9.75, 11.3), each = 4),
  y = rep(4:1, length(table_values)),
  label = unlist(table_values, use.names = FALSE),
  hjust = rep(c(0, rep(0.5, 5)), each = 4),
  fontface = rep(c("bold", rep("plain", 3)), length(table_values))
)
bd_table <- ggplot(table_cells, aes(x = x, y = y, label = label)) +
  geom_text(aes(hjust = hjust, fontface = fontface), size = 3.9) +
  geom_hline(yintercept = c(3.5, 0.5), colour = "#B8C2CC", linewidth = 0.3) +
  scale_x_continuous(limits = c(0, 12), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 4.7), expand = c(0, 0)) +
  labs(
    title = "Boundary density within each illustrative analysis",
    caption = paste(
      "Visium HD: 8-µm bins; Xenium: segmented cells. Expected raw BD uses fixed label counts.",
      "Each row describes one label map; values are not directly comparable across platforms or label sources.",
      sep = "\n"
    )
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.caption = element_text(size = 10, hjust = 0),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
save_plot(
  bd_table,
  file.path(plots_dir, "deOliveira_BD_numeric_summary"),
  7, 2
)

bd$analysis <- paste(bd$platform, bd$label_source, sep = "\n")
# Plot raw observed and expected BD; annotate each panel with relative BD.
bd_long <- rbind(
  data.frame(
    analysis = bd$analysis,
    metric = "Observed",
    value = bd$boundary_density
  ),
  data.frame(
    analysis = bd$analysis,
    metric = "Random-label expectation",
    value = bd$expected_boundary_density
  )
)
bd_long$metric <- factor(
  bd_long$metric, c("Observed", "Random-label expectation")
)
bd_annotations <- data.frame(
  analysis = bd$analysis,
  label = sprintf("Relative BD = %.3f", bd$relative_boundary_density),
  x = 1.5,
  y = 0.99
)

bd_figure <- ggplot(
  bd_long,
  aes(x = metric, y = value, fill = metric)
) +
  geom_col(width = 0.68) +
  geom_text(
    data = bd_annotations,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 1,
    size = 4
  ) +
  facet_wrap(~analysis, nrow = 1) +
  scale_fill_manual(values = c("#3B82F6", "#B8C2CC")) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Boundary density within each illustrative analysis",
    subtitle = paste(
      "Descriptive summaries only; values are not directly comparable",
      "across platforms or label sources"
    ),
    x = NULL, y = "Boundary density", fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )
save_plot(
  bd_figure,
  file.path(plots_dir, "deOliveira_BD_summary"),
  11.5, 5
)

runtime <- rbind(
  read.csv(file.path(results_dir, "deOliveira_VisiumHD_runtime.csv")),
  read.csv(file.path(results_dir, "deOliveira_Xenium_runtime.csv"))
)
# Omit the total row when plotting individual stage timings.
runtime <- runtime[!grepl("^total", runtime$stage), ]
runtime$stage <- gsub("_", " ", runtime$stage, fixed = TRUE)
runtime$stage <- factor(
  runtime$stage,
  levels = rev(unique(runtime$stage))
)
p_runtime <- ggplot(runtime, aes(x = seconds, y = stage, fill = platform)) +
  geom_col(position = "dodge") +
  scale_x_continuous(trans = "log1p") +
  labs(
    title = "Illustrative stage runtimes (single observed run)",
    subtitle = "Descriptive only; platform and location count differ",
    x = "Runtime in seconds (log1p scale)", y = NULL, fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

save_plot(
  p_runtime,
  file.path(plots_dir, "deOliveira_runtime_summary"),
  8, 4.8
)
