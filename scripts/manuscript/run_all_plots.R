################################################################################
# Script to run all plotting scripts in order
################################################################################

# Run from repository root directory.

# GBM
source("scripts/GBM/03_plots/01_plots_GBM_BD.R")
source("scripts/GBM/03_plots/02_plots_GBM_BD_baseline_kmeans.R")
source("scripts/GBM/03_plots/02a_plots_GBM_BD_figure.R")
source("scripts/GBM/03_plots/03_plots_GBM_BD_author_labels.R")
source("scripts/GBM/03_plots/03a_plots_GBM_BD_author_label_clustering.R")
source("scripts/GBM/03_plots/03b_plots_GBM_BD_author_labels_within_tumor.R")
source("scripts/GBM/03_plots/04_plots_GBM_BD_external_clustering.R")
source("scripts/GBM/03_plots/05_plots_GBM_BD_sensitivity_smoothing_method.R")
source("scripts/GBM/03_plots/06_plots_GBM_BD_sensitivity_smoothing_parameter.R")
source("scripts/GBM/03_plots/07_plots_GBM_BD_sensitivity_boundary_k.R")
source("scripts/GBM/03_plots/08_plots_GBM_runtime_comparison.R")

# Kuppe
source("scripts/Kuppe/03_plots/01_plots_Kuppe_BD.R")

# Miyoshi
source("scripts/Miyoshi/03_plots/01_plots_Miyoshi_BD.R")

# Maynard
source("scripts/Maynard/03_plots/01_plots_Maynard_BD.R")

# de Oliveira
source("scripts/deOliveira/03_plots/01_plots_deOliveira.R")

# Runtime
source("scripts/runtime/runtime_summary.R")

# Simulations (run complete scripts)
source("scripts/sim_composition/simulation_composition.R")
source("scripts/sim_oversmoothing/simulation_oversmoothing.R")

# BD neighborhood sensitivity
source("scripts/BD_neighborhood_sensitivity/02_plots_Kuppe_Miyoshi_BD_sensitivity.R")
