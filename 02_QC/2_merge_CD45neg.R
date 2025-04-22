# ---
# author: "Laura Jiménez Gracia"
# date: 2021-29-07
# ---
# This R script merges the cellranger output data [hashed CD45- fraction] into a single Seurat Object.


# Load packages
library(tidyverse)
library(Seurat)

# Paths
path_project_metadata <- here::here("01_cellranger_mapping/data/POST_12_metadata.csv")
path_samples_metadata <- here::here("01_cellranger_mapping/data/POST_12_metadata_demultiplexing.csv")
path_r_objects <- here::here("02_QC/results/R_objects")


# Load metadata information
metadata <- read.csv(path_project_metadata)
metadata_samples <- read.csv(path_samples_metadata)

## Merging metadata_all tables
metadata <- metadata %>% filter(type == "cDNA")
metadata_all <- left_join(metadata, subset(metadata_samples, select = -c(gem_id)), by = "library_name")

## Get list of gem_id and sample_names & merge ids
hashed_gemids <- as.vector(metadata_all$gem_id[metadata_all$hashing == "CellPlex" & metadata_all$fraction == "CD45-"])
sample_ids <- as.vector(metadata_all$sample_id[metadata_all$gem_id %in% hashed_gemids])
hashed_gem_sample_id <- paste0(hashed_gemids, "__", sample_ids)

metadata_all$gem_sample_id <- paste0(metadata_all$gem_id, "__", metadata_all$sample_id)



# Generate a Seurat object
## Get 10X output data paths
path_samplematrix_list <- purrr::map2(hashed_gemids, sample_ids, function(gem_id, sample_id){
  path_sample_matrix <- here::here(paste("01_cellranger_mapping/jobs", gem_id, gem_id, "outs/per_sample_outs", sample_id, "count/sample_feature_bc_matrix", sep = "/"))
  path_sample_matrix
})
path_samplematrix_list <- as.character(path_samplematrix_list)

## Load 10X data and create Seurat object
seurat_obj_hashed_list <- purrr::map2(path_samplematrix_list,
                                      hashed_gem_sample_id,
                                      function(path_sample_matrix, gem_sample_id) {
                                        # Load 10x data
                                        sample_count_matrix <- Read10X(path_sample_matrix)  
                                        # Setup Seurat object, and ignore "Multiplexing Capture" data
                                        seurat_obj <- CreateSeuratObject(counts = sample_count_matrix$`Gene Expression`)
                                        # Add gem_sample_id
                                        seurat_obj$gem_sample_id <- gem_sample_id
                                        seurat_obj
                                      })
names(seurat_obj_hashed_list) <- hashed_gem_sample_id

## Merge Seurat objects
seurat_obj_hashed <- merge(seurat_obj_hashed_list[[1]], 
                           y = seurat_obj_hashed_list[2:length(seurat_obj_hashed_list)],   
                           add.cell.ids = hashed_gemids)
seurat_obj_hashed

## Remove Seurat object lists
rm(seurat_obj_hashed_list)



# Add sample metadata
## Create dataframe with seurat metadata and samples & project
seurat_obj_metadata <- left_join(seurat_obj_hashed@meta.data, metadata_all, by = "gem_sample_id")

## Add metadata information to Seurat object
seurat_obj_hashed$gem_id <- seurat_obj_metadata$gem_id
seurat_obj_hashed$library_name <- seurat_obj_metadata$library_name
seurat_obj_hashed$fraction <- seurat_obj_metadata$fraction
seurat_obj_hashed$sample_id <- seurat_obj_metadata$sample_id
seurat_obj_hashed$genotype <- seurat_obj_metadata$genotype
seurat_obj_hashed$treatment <- seurat_obj_metadata$treatment

## Create new metadata variables by merging
seurat_obj_hashed$condition <- paste(seurat_obj_metadata$genotype, 
                                     seurat_obj_metadata$treatment,
                                     sep = "_")

# Save Seurat object
saveRDS(seurat_obj_hashed, 
        file = paste0(path_r_objects, "/post12_cd45neg_merged.rds"))