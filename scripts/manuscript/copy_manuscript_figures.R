################################################################################
# Script to copy generated figure files to the manuscript directory
################################################################################

# run from repository root directory

figures_dir <- "../smoothclust_paper/figures"
dir.create(figures_dir, showWarnings = FALSE)

# Copy each figure in both PDF and PNG formats when available.
copy_figure <- function(from) {
  stopifnot(file.exists(from))
  extensions <- c("pdf", "png")
  sources <- paste0(tools::file_path_sans_ext(from), ".", extensions)
  available <- file.exists(sources)
  stopifnot(file.copy(
    sources[available], figures_dir, overwrite = TRUE
  ))
}


# GBM --------------------------------------------------------------------------

copy_figure(
  "plots/GBM/baseline_kmeans_comparison/figure/GBM_BD_figure_zh1019t1.png"
)
copy_figure(
  "plots/GBM/external_clustering/figure/GBM_external_clustering_maps_zh1019t1.pdf"
)
copy_figure(
  "plots/GBM/external_clustering/figure/GBM_external_clustering_raw_BD_by_method.pdf"
)
copy_figure(
  "plots/GBM/sensitivity_smoothing_method/figure/GBM_supp_sensitivity_smoothing_method.pdf"
)
copy_figure(
  "plots/GBM/sensitivity_smoothing_parameter/figure/GBM_supp_sensitivity_smoothing_parameter.pdf"
)
copy_figure(
  "plots/GBM/sensitivity_boundary_k/figure/GBM_supp_sensitivity_boundary_k.pdf"
)
copy_figure(
  "plots/GBM/author_labels/supplementary/GBM_author_label_biology_supp.pdf"
)
copy_figure(
  "plots/GBM/author_labels/within_tumor/GBM_author_label_within_tumor_supp.pdf"
)


# Kuppe ------------------------------------------------------------------------

copy_figure(
  "plots/Kuppe/summary_plots/Kuppe_celltype_niche_BD_by_sample_group.pdf"
)
copy_figure(
  "plots/Kuppe/spatial_maps/Kuppe_celltype_niche_example_spatial_maps.pdf"
)
copy_figure(
  "plots/Kuppe/summary_plots/Kuppe_celltype_niche_composition_by_section.pdf"
)
copy_figure(
  "plots/Kuppe/summary_plots/Kuppe_celltype_niche_BD_by_tissue_region.pdf"
)


# Miyoshi ----------------------------------------------------------------------

copy_figure(
  "plots/Miyoshi/summary_plots/Miyoshi_BD_ranked_sections.pdf"
)
copy_figure(
  "plots/Miyoshi/spatial_maps/Miyoshi_cortical_region_BD_all_sections.pdf"
)
# Include the individual PDF/PNG pages of the multipage figure.
for (page in 1:8) {
  copy_figure(sprintf(
    "plots/Miyoshi/spatial_maps/Miyoshi_cortical_region_BD_all_sections_page%02d.pdf",
    page
  ))
}
copy_figure(
  "plots/Miyoshi/spatial_maps/Miyoshi_q9_vs_WM_collapsed_map.pdf"
)
copy_figure(
  "plots/Miyoshi/summary_plots/Miyoshi_BD_q9_vs_WM_collapsed.pdf"
)


# Maynard ----------------------------------------------------------------------

copy_figure(
  "plots/Maynard/summary_plots/Maynard_BD_adjacent_sections.pdf"
)
copy_figure(
  "plots/Maynard/spatial_maps/Maynard_all_sections_manual_layer_maps.pdf"
)
copy_figure(
  "plots/Maynard/spatial_maps/Maynard_all_sections_local_BD_maps.pdf"
)


# de Oliveira ------------------------------------------------------------------

copy_figure(
  "plots/deOliveira/deOliveira_VisiumHD_spatial_summary.pdf"
)
copy_figure(
  "plots/deOliveira/deOliveira_Xenium_spatial_summary.pdf"
)
copy_figure(
  "plots/deOliveira/deOliveira_BD_summary.pdf"
)


# BD calibration simulation ----------------------------------------------------

copy_figure(
  "plots/sim_BD_calibration/sim_BD_calibration.pdf"
)


# Oversmoothing simulation -----------------------------------------------------

copy_figure(
  "plots/sim_oversmoothing/sim_oversmoothing.pdf"
)


# Runtimes ---------------------------------------------------------------------

copy_figure(
  "plots/runtime/runtime_core_workflow_summary.pdf"
)
copy_figure(
  "plots/GBM/runtime/figure/GBM_supp_runtime_summary.pdf"
)
