################################################################################
# Preprocess the Kuppe et al. human myocardial infarction Visium dataset
################################################################################

# Run from the smoothclust-analyses repository root directory.

library(Matrix)
library(readxl)
library(SpatialExperiment)

dataset_id <- "Kuppe"
raw_dir <- file.path("data", dataset_id, "raw")
processed_dir <- file.path("data", dataset_id, "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

clinical <- as.data.frame(read_excel(
  file.path(raw_dir, "Kuppe2022_Supplementary_Table_1.xlsx"), skip = 1
))
names(clinical) <- c(
  "patient_id", "region_novel", "sample_id", "region_code", "age",
  "sex", "sampling_location", "infarction_location",
  "days_after_infarction", "disease", "pathologist_description",
  "donor_information"
)
clinical$days_after_infarction <- suppressWarnings(as.numeric(
  clinical$days_after_infarction
))
clinical <- clinical[!duplicated(clinical$sample_id), ]

visium <- as.data.frame(read_excel(
  file.path(raw_dir, "Kuppe2022_Supplementary_Table_4.xlsx"), skip = 1
))
spaceranger_ids <- aggregate(
  sample_id ~ patient_region_id, visium,
  function(x) paste(x, collapse = "; ")
)
names(spaceranger_ids) <- c("sample_id", "spaceranger_sample_id")

group_labels <- c(
  group_1 = "myogenic-enriched",
  group_2 = "ischaemic-enriched",
  group_3 = "fibrotic-enriched"
)
region_labels <- c(
  CTRL = "Control", RZ = "Remote zone", BZ = "Border zone",
  IZ = "Ischaemic zone", FZ = "Fibrotic zone"
)
prop_names <- c(
  "Adipo", "CM", "Endo", "Fib", "Lymphoid", "Mast", "Myeloid",
  "Neuronal", "PC", "prolif", "vSMCs"
)
prop_columns <- paste0("author_prop_", c(
  "adipocyte", "cardiomyocyte", "endothelial", "fibroblast",
  "lymphoid", "mast", "myeloid", "neuronal", "pericyte",
  "cycling_cells", "vsmc"
))
celltype_names <- c(
  "Adipocyte", "Cardiomyocyte", "Endothelial", "Fibroblast",
  "Lymphoid", "Mast", "Myeloid", "Neuronal", "Pericyte",
  "Cycling cells", "vSMCs"
)

rds_files <- sort(list.files(
  file.path(raw_dir, "hca_rds"), pattern = "[.]rds$", full.names = TRUE
))
count_list <- gene_list <- spot_list <- sample_list <- vector(
  "list", length(rds_files)
)

for (i in seq_along(rds_files)) {
  message("Processing section ", i, "/", length(rds_files), ": ",
          basename(rds_files[i]))
  object <- readRDS(rds_files[i])
  meta <- attr(object, "meta.data")
  assays <- attr(object, "assays")
  counts <- attr(assays[["Spatial"]], "counts")
  props <- attr(assays[["c2l_props"]], "data")
  coords <- attr(attr(object, "images")[[1]], "coordinates")

  barcodes <- rownames(meta)
  stopifnot(
    identical(colnames(counts), barcodes),
    identical(colnames(props), barcodes),
    identical(rownames(coords), barcodes)
  )
  section_id <- unique(as.character(meta$orig.ident))
  sample_id <- unique(as.character(meta$patient_region_id))
  clinical_row <- clinical[match(sample_id, clinical$sample_id), ]

  spot_ids <- paste(section_id, barcodes, sep = "_")
  colnames(counts) <- spot_ids
  count_list[[i]] <- counts
  gene_list[[i]] <- rownames(counts)

  props <- as.matrix(t(props[prop_names, , drop = FALSE]))
  colnames(props) <- prop_columns
  dominant_celltype <- celltype_names[max.col(props, ties.method = "first")]
  author_group_code <- as.character(meta$patient_group)
  region_code <- as.character(meta$major_labl)

  spots <- data.frame(
    spot_id = spot_ids,
    barcode = barcodes,
    dataset_id = dataset_id,
    section_id = section_id,
    sample_id = sample_id,
    patient_id = as.character(meta$patient),
    sample_group = unname(group_labels[author_group_code]),
    author_group_code = author_group_code,
    region_code = region_code,
    region_label = unname(region_labels[region_code]),
    region_novel = clinical_row$region_novel,
    days_after_infarction = clinical_row$days_after_infarction,
    disease = clinical_row$disease,
    age = clinical_row$age,
    sex = clinical_row$sex,
    array_col = coords$col,
    array_row = coords$row,
    author_celltype_niche = as.character(meta$celltype_niche),
    author_molecular_niche = as.character(meta$molecular_niche),
    author_dominant_celltype = dominant_celltype,
    author_n_counts = meta$nCount_Spatial,
    author_n_genes = meta$nFeature_Spatial,
    author_percent_mito = meta$percent.mt,
    stringsAsFactors = FALSE
  )
  spots[prop_columns] <- props
  spot_list[[i]] <- spots

  sample_list[[i]] <- data.frame(
    section_id = section_id,
    sample_id = sample_id,
    patient_id = unique(as.character(meta$patient)),
    sample_group = unique(unname(group_labels[author_group_code])),
    author_group_code = unique(author_group_code),
    region_code = unique(region_code),
    region_label = unique(unname(region_labels[region_code])),
    region_novel = clinical_row$region_novel,
    days_after_infarction = clinical_row$days_after_infarction,
    disease = clinical_row$disease,
    age = clinical_row$age,
    sex = clinical_row$sex,
    sampling_location = clinical_row$sampling_location,
    infarction_location = clinical_row$infarction_location,
    pathologist_description = clinical_row$pathologist_description,
    donor_information = clinical_row$donor_information,
    spaceranger_sample_id = spaceranger_ids$spaceranger_sample_id[
      match(sample_id, spaceranger_ids$sample_id)
    ],
    hca_rds_file = basename(rds_files[i]),
    n_spots = nrow(meta),
    n_celltype_niches = length(unique(meta$celltype_niche)),
    n_molecular_niches = length(unique(meta$molecular_niche)),
    stringsAsFactors = FALSE
  )
  rm(object, assays, counts, props)
  gc(FALSE)
}

genes <- Reduce(intersect, gene_list)
message("Retaining ", length(genes), " genes shared across all sections")
count_list <- Map(
  function(x, g) x[match(genes, g), , drop = FALSE], count_list, gene_list
)
counts <- do.call(cbind, count_list)
spots <- do.call(rbind, spot_list)
sample_data <- do.call(rbind, sample_list)
rownames(spots) <- spots$spot_id

spatial_coords <- cbind(
  x = spots$array_col / 2,
  y = spots$array_row * sqrt(3) / 2
)
rownames(spatial_coords) <- spots$spot_id
gene_data <- DataFrame(
  gene_name = genes,
  symbol = genes,
  is_mito = grepl("^MT-", genes),
  row.names = genes
)
spe <- SpatialExperiment(
  assays = list(counts = counts),
  rowData = gene_data,
  colData = DataFrame(spots),
  spatialCoords = spatial_coords
)

metadata(spe)$dataset_id <- dataset_id
metadata(spe)$hca_project <- "e9f36305-d857-44a3-93f0-df4e6007dc97"
metadata(spe)$source <- paste(
  "Kuppe et al. (2022), doi:10.1038/s41586-022-05060-x"
)
metadata(spe)$annotation_description <- paste(
  "author_celltype_niche and author_molecular_niche are exact deposited",
  "spot assignments; author_prop_* are deposited cell2location proportions"
)
metadata(spe)$counts_description <- paste(
  "raw counts restricted to genes shared across all 28 Visium sections"
)
metadata(spe)$coordinate_description <- paste(
  "author Visium array positions transformed to an equilateral hexagonal",
  "lattice"
)
metadata(spe)$section_description <- paste(
  "28 Visium sections from 27 patient-region samples; RZ_P6 has two sections"
)
metadata(spe)$sample_metadata <- sample_data

saveRDS(spe, file.path(processed_dir, "spe_Kuppe.rds"))
write.csv(
  sample_data, file.path(processed_dir, "Kuppe_sample_metadata.csv"),
  row.names = FALSE
)
