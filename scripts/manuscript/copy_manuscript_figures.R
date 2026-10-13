################################################################################
# Script to copy generated figure files to the manuscript directory
################################################################################

# run from repository root directory

# Current final figures, exported in both PDF and PNG formats.
figures <- c(
  # GBM
  "plots/GBM/baseline_kmeans_comparison/figure/GBM_BD_figure_zh1019t1.pdf",
  "plots/GBM/external_clustering/figure/GBM_external_clustering_maps_zh1019t1.pdf",
  "plots/GBM/external_clustering/figure/GBM_external_clustering_raw_BD_by_method.pdf",
  "plots/GBM/sensitivity_smoothing_method/figure/GBM_supp_sensitivity_smoothing_method.pdf",
  "plots/GBM/sensitivity_smoothing_parameter/figure/GBM_supp_sensitivity_smoothing_parameter.pdf",
  "plots/GBM/sensitivity_boundary_k/figure/GBM_supp_sensitivity_boundary_k.pdf",
  "plots/GBM/author_labels/supplementary/GBM_author_label_biology_supp.pdf",
  "plots/GBM/author_labels/within_tumor/GBM_author_label_within_tumor_supp.pdf",
  # Kuppe
  "plots/Kuppe/summary_plots/Kuppe_celltype_niche_BD_by_sample_group.pdf",
  "plots/Kuppe/spatial_maps/Kuppe_celltype_niche_example_spatial_maps.pdf",
  "plots/Kuppe/summary_plots/Kuppe_celltype_niche_composition_by_section.pdf",
  # Miyoshi
  "plots/Miyoshi/summary_plots/Miyoshi_BD_ranked_sections.pdf",
  "plots/Miyoshi/spatial_maps/Miyoshi_q9_vs_WM_collapsed_map.pdf",
  "plots/Miyoshi/summary_plots/Miyoshi_BD_q9_vs_WM_collapsed.pdf",
  # Maynard
  "plots/Maynard/spatial_maps/Maynard_all_sections_manual_layer_maps.pdf",
  "plots/Maynard/spatial_maps/Maynard_all_sections_local_BD_maps.pdf",
  # de Oliveira
  "plots/deOliveira/deOliveira_VisiumHD_spatial_summary.pdf",
  "plots/deOliveira/deOliveira_Xenium_spatial_summary.pdf",
  "plots/deOliveira/deOliveira_BD_numeric_summary.pdf",
  # Simulations
  "plots/sim_composition/sim_composition.pdf",
  "plots/sim_oversmoothing/sim_oversmoothing.pdf",
  # Runtime
  "plots/runtime/runtime_core_workflow_summary.pdf"
)

figures <- c(
  figures,
  sub("\\.pdf$", ".png", figures),
  # Complete multipage atlas, exported only as PDF.
  "plots/Miyoshi/spatial_maps/Miyoshi_cortical_region_BD_all_sections.pdf"
)

missing_figures <- figures[!file.exists(figures)]
if (length(missing_figures)) {
  stop("Missing required final figures:\n",
       paste(missing_figures, collapse = "\n"), call. = FALSE)
}

figures_dir <- "../smoothclust_paper/figures"
dir.create(figures_dir, showWarnings = FALSE)

stopifnot(file.copy(figures, figures_dir, overwrite = TRUE))
