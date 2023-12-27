# set up object and calculate HVGs

library(STexampleData)
spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]

library(scran)
spe <- logNormCounts(spe)

is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
spe <- spe[!is_mito, ]

dec <- modelGeneVar(spe)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
spe_hvgs <- spe[top_hvgs, ]

dim(spe_hvgs)


# calculate neighbors (including self)

colData(spe_hvgs)$neighbors <- as.list(rep(NA, ncol(spe_hvgs)))

arrayrow <- colData(spe_hvgs)$array_row
arraycol <- colData(spe_hvgs)$array_col

# slow step (can improve)
for (i in 1:ncol(spe_hvgs)) {
  neighbors <- i
  for (j in setdiff(1:ncol(spe_hvgs), i)) {
    # defining neighbors as within row and column distance threshold
    if (abs(arrayrow[i] - arrayrow[j]) <= 5 & (abs(arraycol[i] - arraycol[j])) <= 5) {
      neighbors <- c(neighbors, j)
    }
  }
  # neighbors stored as indices
  colData(spe_hvgs)$neighbors[[i]] <- neighbors
  print(i)
}

head(colData(spe_hvgs))
head(colData(spe_hvgs)$neighbors)

# check
colData(spe_hvgs)$neighbors[[1]]


# calculate average logcounts across neighbors (vectorized for faster runtime)
# (alternatively: use kernels for more sophisticated approach)

smoothed_logcounts <- matrix(NA, nrow = nrow(spe_hvgs), ncol = ncol(spe_hvgs))

for (i in 1:ncol(smoothed_logcounts)) {
  smoothed_logcounts[, i] <- rowMeans(logcounts(spe_hvgs)[, colData(spe_hvgs)$neighbors[[i]], drop = FALSE])
  print(i)
}

rownames(smoothed_logcounts) <- rownames(spe_hvgs)
colnames(smoothed_logcounts) <- colnames(spe_hvgs)

assays(spe_hvgs)[["smoothed_logcounts"]] <- smoothed_logcounts

assayNames(spe_hvgs)


### too slow (nested loop)

# smoothed_logcounts <- matrix(NA, nrow = nrow(spe_hvgs), ncol = ncol(spe_hvgs))
# 
# for (g in 1:nrow(smoothed_logcounts)) {
#   for (j in 1:ncol(smoothed_logcounts)) {
#     smoothed_logcounts[g, j] <- mean(logcounts(spe_hvgs)[g, colData(spe_hvgs)$neighbors[[j]]])
#     print(j)
#   }
#   print(g)
# }


# plots for checking (PCP4 gene)

library(ggplot2)

df <- cbind(as.data.frame(spatialCoords(spe_hvgs)), 
            gene = smoothed_logcounts["ENSG00000183036", ], 
            gene_original = logcounts(spe)[33335, ])

colnames(df) <- c("x", "y", "gene", "gene_original")

ggplot(df, aes(x = x, y = y, color = gene)) + geom_point(size = 0.5) + scale_y_reverse() + coord_fixed()
ggplot(df, aes(x = x, y = y, color = gene_original)) + geom_point(size = 0.5) + scale_y_reverse() + coord_fixed()


# dimensionality reduction and clustering

library(scater)
library(scran)

# compute PCA
set.seed(123)
spe_hvgs <- runPCA(spe_hvgs, exprs_values = "smoothed_logcounts")
# compute UMAP on top 50 PCs
set.seed(123)
spe_hvgs <- runUMAP(spe_hvgs, dimred = "PCA")
colnames(reducedDim(spe_hvgs, "UMAP")) <- paste0("UMAP", 1:2)


# graph-based clustering
set.seed(123)
#k <- 10
k <- 100
g <- buildSNNGraph(spe_hvgs, k = k, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
table(clus)
colLabels(spe_hvgs) <- factor(clus)


# alternatively: k-means clustering (works better)
k <- 6
set.seed(100)
clust <- kmeans(reducedDim(spe_hvgs, "PCA"), centers = k)$cluster
table(clust)
colLabels(spe_hvgs) <- factor(clust)


# plots
library(ggspavis)

plotSpots(spe_hvgs, annotate = "label", palette = "libd_layer_colors")
plotSpots(spe_hvgs, annotate = "label", palette = unname(palette.colors(36, "Polychrome 36")))

plotSpots(spe_hvgs, annotate = "ground_truth", palette = "libd_layer_colors")


# ARI

library(mclust)
adjustedRandIndex(colData(spe_hvgs)$label, colData(spe)$ground_truth)

