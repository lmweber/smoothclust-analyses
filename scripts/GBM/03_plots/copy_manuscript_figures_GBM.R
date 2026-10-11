################################################################################
# Script to copy generated figure files to the manuscript repository
################################################################################

# run from repository root directory

stopifnot(
  file.copy("plots/GBM/baseline_kmeans_comparison/figure/GBM_BD_figure_zh1019t1.png", 
            "../smoothclust_paper/figures/GBM_BD_figure_zh1019t1.png", 
            overwrite = TRUE))
stopifnot(
  file.copy("plots/GBM/external_clustering/figure/GBM_external_clustering_raw_BD_by_method.pdf", 
            "../smoothclust_paper/figures/GBM_external_clustering_raw_BD_by_method.pdf", 
            overwrite = TRUE))
stopifnot(
  file.copy("plots/GBM/external_clustering/figure/GBM_external_clustering_maps_zh1019t1.pdf", 
            "../smoothclust_paper/figures/GBM_external_clustering_maps_zh1019t1.pdf", 
            overwrite = TRUE))
stopifnot(
  file.copy("plots/GBM/sensitivity_smoothing_method/figure/GBM_supp_sensitivity_smoothing_method.pdf", 
            "../smoothclust_paper/figures/GBM_supp_sensitivity_smoothing_method.pdf", 
            overwrite = TRUE))
stopifnot(
  file.copy("plots/GBM/sensitivity_smoothing_parameter/figure/GBM_supp_sensitivity_smoothing_parameter.pdf", 
            "../smoothclust_paper/figures/GBM_supp_sensitivity_smoothing_parameter.pdf", 
            overwrite = TRUE))
stopifnot(
  file.copy("plots/GBM/sensitivity_boundary_k/figure/GBM_supp_sensitivity_boundary_k.pdf", 
            "../smoothclust_paper/figures/GBM_supp_sensitivity_boundary_k.pdf", 
            overwrite = TRUE))
stopifnot(
  file.copy("plots/GBM/author_labels/supplementary/GBM_author_label_biology_supp.pdf", 
            "../smoothclust_paper/figures/GBM_author_label_biology_supp.pdf", 
            overwrite = TRUE))
stopifnot(
  file.copy("plots/GBM/author_labels/within_tumor/GBM_author_label_within_tumor_supp.pdf",
            "../smoothclust_paper/figures/GBM_author_label_within_tumor_supp.pdf",
            overwrite = TRUE))
stopifnot(
  file.copy("plots/GBM/runtime/figure/GBM_supp_runtime_summary.pdf", 
            "../smoothclust_paper/figures/GBM_supp_runtime_summary.pdf", 
            overwrite = TRUE))
