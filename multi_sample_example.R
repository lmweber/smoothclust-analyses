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


# simple normalization by sample (by mean sum UMIs across spots per sample)

library(tidyverse)

col_sums <- colSums(counts(spe))
sample_ids <- colData(spe)$sample_id

df <- data.frame(col_sums, sample_ids) |> 
  group_by(sample_ids) |> 
  summarize(mean_by_sample = mean(col_sums))

scale_factors_by_sample <- df$mean_by_sample / mean(df$mean_by_sample)
scale_factors_by_sample

n_spots_per_sample <- unname(table(colData(spe)$sample_id))
scale_factors_by_sample_vec <- rep(scale_factors_by_sample, times = n_spots_per_sample)

# scale counts
# skip this if not scaling (e.g. using Harmony later)
counts(spe) <- t(t(counts(spe)) / scale_factors_by_sample_vec)


# run smoothclust

library(smoothclust)

# run on raw counts and re-calculate logcounts later
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
assayNames(spe)


# run smoothclust individually per sample
for (i in seq_along(sample_names)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_names[i]]
  
  print(sample_names[i])
  print(dim(spe_sub))
  
  print(Sys.time())
  spe_sub <- smoothclust(spe_sub, assay_name = "counts", 
                         method = "uniform", bandwidth = 0.05)
  print(Sys.time())
  
  # remove normcounts assay to reduce memory usage
  print(assayNames(spe_sub))
  assay(spe_sub, "counts") <- NULL
  print(assayNames(spe_sub))
  
  spe_list[[i]] <- spe_sub
}


# stack into combined SPE object
# note: large amount of memory needed

# slow (~30 sec)
Sys.time()
spe_combined <- do.call("cbind", spe_list)
Sys.time()

dim(spe_combined)
assayNames(spe_combined)
table(colData(spe_combined)$sample_id)


# calculate logcounts
# note: spot-level normalization, in addition to any previous sample-level normalization
library(scater)
library(scran)
spe_combined <- logNormCounts(spe_combined, assay.type = "counts_smooth")
assayNames(spe_combined)

# remove parts of object to reduce memory usage
assay(spe_combined, "counts_smooth") <- NULL
# remove objects to reduce memory usage
rm(spe)
rm(spe_list)
rm(spe_sub)


# clustering

# preprocessing steps for clustering
# remove mitochondrial genes
is_mito <- grepl("(^mt-)", rowData(spe_combined)$gene_name, ignore.case = TRUE)
table(is_mito)
# slow (~10 sec)
Sys.time()
spe_combined <- spe_combined[!is_mito, ]
Sys.time()
dim(spe_combined)

# select top highly variable genes (HVGs)
# slow / high memory usage
# need to run on compute cluster instead due to memory requirements
# use previous set of HVGs for now
# library(BiocParallel)
# Sys.time()
# dec <- modelGeneVar(spe_combined, 
#                     block = as.factor(colData(spe_combined)$sample_id), 
#                     BPPARAM = MulticoreParam(workers = 6))
# Sys.time()
# top_hvgs <- getTopHVGs(dec, prop = 0.1)

# re-use previous set of HVGs for now
table(rowData(spe_combined)$is_top_hvg)
top_hvgs <- rowData(spe_combined)[rowData(spe_combined)$is_top_hvg, "gene_id"]
length(top_hvgs)
spe_combined <- spe_combined[top_hvgs, ]
dim(spe_combined)

# dimensionality reduction
# compute PCA on top HVGs
# slow (~1-2 min)
set.seed(123)
Sys.time()
spe_combined <- runPCA(spe_combined)
Sys.time()

# sample-level normalization using Harmony
library(harmony)
meta_data <- scale_factors_by_sample_vec
names(meta_data) <- sample_ids
harmony_embeddings <- RunHarmony(reducedDim(spe_combined, "PCA"), meta_data, "dataset")

reducedDim(spe_combined, "HARM") <- harmony_embeddings


# run k-means clustering (for selected number of clusters)
set.seed(123)
k <- 2
#clus <- kmeans(reducedDim(spe_combined, "PCA"), centers = k)$cluster
clus <- kmeans(reducedDim(spe_combined, "HARM"), centers = k)$cluster
table(clus)
colLabels(spe_combined) <- factor(clus)

# alternatively: run graph-based clustering
# set.seed(123)
# k <- 20
# #g <- buildSNNGraph(spe_combined, k = k, use.dimred = "PCA")
# g <- buildSNNGraph(spe_combined, k = k, use.dimred = "HARM")
# g_walk <- igraph::cluster_walktrap(g)
# clus <- g_walk$membership
# table(clus)
# colLabels(spe_combined) <- factor(clus)


# plot clustering

# segments white matter quite well when using 2 clusters

# clusters split by sample when using more clusters
# using Harmony embeddings merges clusters across samples quite well
# note: may need slightly larger bandwidth for smoother boundaries

# note: re-calculate HVGs using higher memory on compute cluster
# or calculate them individually per sample and then merge the vectors

plotVisium(spe_combined, annotate = "label", 
           facets = "sample_id")#, pal = "libd_layer_colors")

plotVisium(spe_combined, annotate = "label", 
           facets = "sample_id", pal = "libd_layer_colors")

plotVisium(spe_combined, annotate = "layer_guess_reordered_short", 
           facets = "sample_id", pal = "libd_layer_colors")
