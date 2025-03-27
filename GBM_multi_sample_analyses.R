# Multi-sample example using GBM dataset

library(SpatialExperiment)
library(ggspavis)
library(smoothclust)
library(scater)
library(scran)


# load data

dir <- "../data/GBM/GBM_ZH881T1"
fnm <- file.path(dir, "raw_feature_bc_matrix")
sce <- DropletUtils::read10xCounts(fnm)
img <- readImgData(
  path = file.path(dir, "spatial"), 
  sample_id = "ZH881T1")
fnm <- file.path(dir, "spatial", "tissue_positions_list.csv")
xyz <- read.csv(fnm, header = FALSE, 
                col.names = c(
                  "barcode", "in_tissue", "array_row", "array_col", 
                  "pxl_row_in_fullres", "pxl_col_in_fullres"))
xyz <- xyz[xyz$in_tissue == 1, ]
rd <- S4Vectors::DataFrame(
  symbol = rowData(sce)$Symbol)
spe_ZH881T1 <- SpatialExperiment(
  assays = list(counts = assay(sce)), 
  rowData = rd, 
  colData = DataFrame(xyz), 
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"), 
  imgData = img, 
  sample_id = "ZH881T1")


dir <- "../data/GBM/GBM_ZH8812bulk"
fnm <- file.path(dir, "raw_feature_bc_matrix")
sce <- DropletUtils::read10xCounts(fnm)
img <- readImgData(
  path = file.path(dir, "spatial"), 
  sample_id = "ZH8812bulk")
fnm <- file.path(dir, "spatial", "tissue_positions_list.csv")
xyz <- read.csv(fnm, header = FALSE, 
                col.names = c(
                  "barcode", "in_tissue", "array_row", "array_col", 
                  "pxl_row_in_fullres", "pxl_col_in_fullres"))
xyz <- xyz[xyz$in_tissue == 1, ]
rd <- S4Vectors::DataFrame(
  symbol = rowData(sce)$Symbol)
spe_ZH8812bulk <- SpatialExperiment(
  assays = list(counts = assay(sce)), 
  rowData = rd, 
  colData = DataFrame(xyz), 
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"), 
  imgData = img, 
  sample_id = "ZH8812bulk")


dir <- "../data/GBM/GBM_ZH8811Abulk"
fnm <- file.path(dir, "raw_feature_bc_matrix")
sce <- DropletUtils::read10xCounts(fnm)
img <- readImgData(
  path = file.path(dir, "spatial"), 
  sample_id = "ZH8811Abulk")
fnm <- file.path(dir, "spatial", "tissue_positions_list.csv")
xyz <- read.csv(fnm, header = FALSE, 
                col.names = c(
                  "barcode", "in_tissue", "array_row", "array_col", 
                  "pxl_row_in_fullres", "pxl_col_in_fullres"))
xyz <- xyz[xyz$in_tissue == 1, ]
rd <- S4Vectors::DataFrame(
  symbol = rowData(sce)$Symbol)
spe_ZH8811Abulk <- SpatialExperiment(
  assays = list(counts = assay(sce)), 
  rowData = rd, 
  colData = DataFrame(xyz), 
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"), 
  imgData = img, 
  sample_id = "ZH8811Abulk")


# combine into single object

dim(spe_ZH881T1)
dim(spe_ZH8812bulk)
dim(spe_ZH8811Abulk)

colData(spe_ZH881T1)
colData(spe_ZH8812bulk)
colData(spe_ZH8811Abulk)

rowData(spe_ZH881T1)
rowData(spe_ZH8812bulk)
rowData(spe_ZH8811Abulk)


spe <- cbind(spe_ZH881T1, spe_ZH8812bulk, spe_ZH8811Abulk)

dim(spe)
colData(spe)
rowData(spe)

colData(spe)$sample_id <- factor(colData(spe)$sample_id, 
                                 levels = c("ZH881T1", "ZH8812bulk", "ZH8811Abulk"))

dim(spe)
colData(spe)
rowData(spe)
table(colData(spe)$sample_id)

# plot points by sample
plotVisium(spe)


# check distributions of UMI counts per sample

library(tidyverse)

colData(spe)$sum_umi <- colSums(counts(spe))

df <- as.data.frame(colData(spe))
head(df)

# plot by sample
ggplot(df, aes(x = sample_id, y = sum_umi, color = sample_id)) + 
  geom_jitter(alpha = 0.1, width = 0.25) + 
  geom_boxplot(color = "black", fill = NA, width = 0.5, outlier.size = 0.1) + 
  theme_bw()

ggsave("plots/GBM_UMI_by_sample.png", width = 5, height = 4.5)


# remove genes with low expression (for smoothclust)

dim(spe)
is_low <- rowSums(counts(spe)) <= 100
table(is_low)

spe <- spe[!is_low, ]

dim(spe)
colData(spe)
rowData(spe)


# run smoothclust

# run on raw counts and re-calculate logcounts later
assayNames(spe)

sample_names <- names(table(colData(spe)$sample_id))
sample_names

spe_list <- as.list(rep(NA, length(sample_names)))
names(spe_list) <- sample_names
spe_list
length(spe_list)

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
# note: large amount of memory needed - to do: improve this

# slow (~ up to 30 sec)
Sys.time()
spe_combined <- do.call("cbind", spe_list)
Sys.time()

dim(spe_combined)
assayNames(spe_combined)
colData(spe_combined)
rowData(spe_combined)
table(colData(spe_combined)$sample_id)

# calculate logcounts
spe_combined <- logNormCounts(spe_combined, assay.type = "counts_smooth")
assayNames(spe_combined)

dim(spe_combined)
assayNames(spe_combined)
colData(spe_combined)
rowData(spe_combined)


# clustering

# preprocessing steps for clustering
# remove mitochondrial genes
is_mito <- grepl("(^mt-)", rowData(spe_combined)$symbol, ignore.case = TRUE)
table(is_mito)
# slow (~ up to 10 sec)
Sys.time()
spe_combined <- spe_combined[!is_mito, ]
Sys.time()
dim(spe_combined)

# select top highly variable genes (HVGs)
# note: slow / high memory usage (~30 sec)
library(BiocParallel)
Sys.time()
dec <- modelGeneVar(spe_combined,
                    block = as.factor(colData(spe_combined)$sample_id), 
                    BPPARAM = MulticoreParam(workers = 14))
Sys.time()
top_hvgs <- getTopHVGs(dec, prop = 0.1)

length(top_hvgs)
head(top_hvgs)


# dimensionality reduction
# compute PCA on top HVGs
# slow (~ up to 1-2 min)
set.seed(123)
Sys.time()
spe_combined <- runPCA(spe_combined)
Sys.time()

# sample-level normalization using Harmony
library(harmony)
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
head(scale_factors_by_sample_vec)
length(scale_factors_by_sample_vec)
dim(spe)
meta_data <- scale_factors_by_sample_vec
names(meta_data) <- sample_ids

harmony_embeddings <- RunHarmony(reducedDim(spe_combined, "PCA"), meta_data, "dataset")

reducedDim(spe_combined, "HARM") <- harmony_embeddings


# run k-means clustering (for selected number of clusters)
set.seed(123)
k <- 3
# use either PCA or HARM embeddings as input for clustering
# clus <- kmeans(reducedDim(spe_combined, "PCA"), centers = k)$cluster
clus <- kmeans(reducedDim(spe_combined, "HARM"), centers = k)$cluster
table(clus)
colLabels(spe_combined) <- factor(clus)


# plot clustering

plotVisium(spe_combined, annotate = "label", 
           facets = "sample_id", pal = c("goldenrod1", "cornflowerblue", "tomato"))

plotVisium(spe_combined, annotate = "label", 
           facets = "sample_id", pal = "libd_layer_colors")

ggsave("plots/smoothclust_GBM_3samples.png", width = 6, height = 6)


# libd_layer_colors <-  c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", 
#                         "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
