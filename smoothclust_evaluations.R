# Evaluations for smoothclust algorithm
# Lukas Weber, Dec 2023

library(smoothclust)
library(STexampleData)
library(scran)
library(scater)
library(ggplot2)
library(ggspavis)
library(mclust)


# set up object and calculate HVGs

spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]

spe <- logNormCounts(spe)

is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
spe <- spe[!is_mito, ]

# keep full object for plotting
spe_full <- spe

dec <- modelGeneVar(spe)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
spe <- spe[top_hvgs, ]

dim(spe)


# run smoothclust
# runtime: ~1 min with default parameters

spe <- smoothclust(spe)


# plots for checking (PCP4 gene)

df <- cbind(as.data.frame(spatialCoords(spe)), 
            gene = assay(spe, "logcounts_smooth")["ENSG00000183036", ], 
            gene_original = logcounts(spe_full)[33335, ])

colnames(df) <- c("x", "y", "gene", "gene_original")

ggplot(df, aes(x = x, y = y, color = gene)) + geom_point(size = 1.25) + scale_y_reverse() + coord_fixed()
ggplot(df, aes(x = x, y = y, color = gene_original)) + geom_point(size = 1.25) + scale_y_reverse() + coord_fixed()


# dimensionality reduction and clustering

# compute PCA
set.seed(123)
spe <- runPCA(spe, exprs_values = "logcounts_smooth")
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# # graph-based clustering
# set.seed(123)
# #k <- 10
# k <- 100
# g <- buildSNNGraph(spe, k = k, use.dimred = "PCA")
# g_walk <- igraph::cluster_walktrap(g)
# clus <- g_walk$membership
# table(clus)
# colLabels(spe) <- factor(clus)

# alternatively: k-means clustering (works better)
# (note 'logcounts_smooth' assay is not sparse, vs. 'logcounts' which is sparse)
k <- 6
set.seed(100)
clust <- kmeans(reducedDim(spe, "PCA"), centers = k)$cluster
table(clust)
colLabels(spe) <- factor(clust)


# plots

plotSpots(spe, annotate = "label", palette = "libd_layer_colors")
plotSpots(spe, annotate = "label", palette = unname(palette.colors(36, "Polychrome 36")))

plotSpots(spe, annotate = "ground_truth", palette = "libd_layer_colors")


# note bandwidth = 0.05 is slightly too large for this dataset since layer 1 disappears


# ARI

adjustedRandIndex(colData(spe)$label, colData(spe_full)$ground_truth)

