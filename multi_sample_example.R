# Clustering example for multi-sample data
# Lukas Weber, 2024-03-30

# run on compute cluster due to memory usage


# download data

library(spatialLIBD)
spe <- fetch_data(type = "spe")

dim(spe)
head(colData(spe), 2)
table(colData(spe)$sample_id)
table(colData(spe)$layer_guess_reordered_short)


# plot

library(ggspavis)

plotVisium(spe, annotate = "layer_guess_reordered_short", 
           facets = "sample_id", pal = "libd_layer_colors") + 
  guides(fill = guide_legend(title = "layer"))


# check distributions of UMI counts per sample

library(tidyverse)

df <- as.data.frame(colData(spe)) |> 
  select(c("sample_id", "layer_guess_reordered_short", "sum_umi")) |> 
  mutate(layer = layer_guess_reordered_short)

head(df)

# plot by sample
ggplot(df, aes(x = sample_id, y = sum_umi, color = sample_id)) + 
  geom_jitter(alpha = 0.1) + 
  geom_boxplot(color = "black", fill = NA) + 
  theme_bw()

# plot by sample and layer
ggplot(df, aes(x = sample_id, y = sum_umi, color = layer)) + 
  geom_boxplot() + 
  theme_bw()


# run smoothclust
# note: no normalization across samples

library(smoothclust)

# run on raw counts and re-calculate logcounts
assayNames(spe)

# run separately per sample
sample_names <- names(table(colData(spe)$sample_id))
sample_names

spe_list <- as.list(rep(NA, length(sample_names)))
names(spe_list) <- sample_names
spe_list

# remove parts of object to reduce memory usage on laptop
logcounts(spe) <- NULL
reducedDims(spe) <- NULL

# assayNames(spe)
# dim(spe)
# head(colData(spe), 2)
# head(rowData(spe), 2)
# 
# spe <- spe[rowData(spe)$is_top_hvg, ]
# dim(spe)
# assay(spe, "logcounts") <- NULL
# colData(spe) <- colData(spe)[, c("sample_id", "sum_umi", "layer_guess_reordered_short")]
# rowData(spe) <- rowData(spe)[, c("gene_id", "gene_name")]
# reducedDims(spe) <- NULL
# 
# assayNames(spe)
# head(colData(spe), 2)
# head(rowData(spe), 2)
# dim(spe)
# 
# counts(spe)[1:6, 1:6]


# run smoothclust individually per sample
for (i in seq_along(sample_names)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_names[i]]
  
  print(sample_names[i])
  print(dim(spe_sub))
  
  print(Sys.time())
  spe_sub <- smoothclust(spe_sub, assay_name = "counts", 
                         method = "uniform", bandwidth = 0.05)
  print(Sys.time())
  
  # remove counts assay to reduce memory usage
  print(assayNames(spe_sub))
  assay(spe_sub, "counts") <- NULL
  print(assayNames(spe_sub))
  
  spe_list[[i]] <- spe_sub
}


# stack into combined SPE object
# note: large amount of memory needed

spe_combined <- do.call("cbind", spe_list)

dim(spe_combined)
assayNames(spe_combined)
table(colData(spe_combined)$sample_id)


# re-calculate logcounts
library(scater)
library(scran)
spe_combined <- logNormCounts(spe_combined, assay.type = "counts_smooth")
assayNames(spe_combined)


# clustering

# preprocessing steps for clustering
# remove mitochondrial genes
is_mito <- grepl("(^mt-)", rowData(spe_combined)$gene_name, ignore.case = TRUE)
table(is_mito)
spe_combined <- spe_combined[!is_mito, ]
dim(spe_combined)

# select top highly variable genes (HVGs)
dec <- modelGeneVar(spe_combined)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
spe_combined <- spe_combined[top_hvgs, ]
dim(spe_combined)

# dimensionality reduction
# compute PCA on top HVGs
set.seed(123)
spe_combined <- runPCA(spe_combined)

# run k-means clustering (for selected number of clusters)
set.seed(123)
k <- 2
clus <- kmeans(reducedDim(spe_combined, "PCA"), centers = k)$cluster
table(clus)
colLabels(spe_combined) <- factor(clus)

# alternatively: run graph-based clustering
# set.seed(123)
# k <- 20
# g <- buildSNNGraph(spe_combined, k = k, use.dimred = "PCA")
# g_walk <- igraph::cluster_walktrap(g)
# clus <- g_walk$membership
# table(clus)
# colLabels(spe_combined) <- factor(clus)


# plot clustering

plotVisium(spe_combined, annotate = "label", 
           facets = "sample_id")#, pal = "libd_layer_colors")

plotVisium(spe_combined, annotate = "layer_guess_reordered_short", 
           facets = "sample_id", pal = "libd_layer_colors")