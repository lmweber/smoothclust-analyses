# Single-sample example using GBM dataset

library(SpatialExperiment)
library(ggspavis)
library(smoothclust)
library(scater)
library(scran)


# load data
# select one of these three samples

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

spe <- spe_ZH881T1
spe <- spe_ZH8812bulk
spe <- spe_ZH8811Abulk

dim(spe)
colData(spe)
rowData(spe)


# plot points
plotVisium(spe)


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

spe <- smoothclust(spe, assay_name = "counts", 
                   method = "uniform", bandwidth = 0.05)


# calculate logcounts
spe <- logNormCounts(spe, assay.type = "counts_smooth")
assayNames(spe)


# clustering

# preprocessing steps for clustering
# remove mitochondrial genes
is_mito <- grepl("(^mt-)", rowData(spe)$symbol, ignore.case = TRUE)
table(is_mito)
# slow (~ up to 10 sec)
Sys.time()
spe <- spe[!is_mito, ]
Sys.time()
dim(spe)

# select top highly variable genes (HVGs)
Sys.time()
dec <- modelGeneVar(spe)
Sys.time()
top_hvgs <- getTopHVGs(dec, prop = 0.1)

length(top_hvgs)
head(top_hvgs)


# dimensionality reduction
# compute PCA on top HVGs
# slow (~ up to 1-2 min)
set.seed(123)
Sys.time()
spe <- runPCA(spe)
Sys.time()

# run k-means clustering (for selected number of clusters)
set.seed(123)
k <- 3
# use PCA (single sample)
clus <- kmeans(reducedDim(spe, "PCA"), centers = k)$cluster
table(clus)
colLabels(spe) <- factor(clus)


# plot clustering

plotVisium(spe, annotate = "label", 
           facets = "sample_id", pal = c("goldenrod1", "cornflowerblue", "tomato"))

plotVisium(spe, annotate = "label", 
           facets = "sample_id", pal = "libd_layer_colors")

ggsave("plots/smoothclust_GBM_1sample_ZH881T1.png", width = 6, height = 6)
ggsave("plots/smoothclust_GBM_1sample_ZH8812bulk.png", width = 6, height = 6)
ggsave("plots/smoothclust_GBM_1sample_ZH8811Abulk.png", width = 6, height = 6)


# libd_layer_colors <-  c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", 
#                         "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
