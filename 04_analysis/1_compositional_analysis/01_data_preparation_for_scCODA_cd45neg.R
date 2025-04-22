# ---
# author: "Laura Jiménez Gracia"
# date: 2025-03-03
# This R script allows to prepare 


# Pre-processing
## Load packages
library(tidyverse)
library(Seurat)

## Paths
path_r_tables_in <- here::here("03_clustering_annotation/results/tables")
path_r_tables_out <- here::here("04_analysis/1_compositional_analysis/results/tables")

## Load data
seurat_obj_metadata <- read.csv(paste0(path_r_tables_in, "/post12_cd45neg_alldata_cleaned_final_metadata.csv"))
seurat_obj_metadata

# Compute number of cells 

## Cell-lineage & sample (treatment)
celllineage_prop_df <- seurat_obj_metadata %>% 
  select(c("sample_id", "celltypes_1")) %>%
  dplyr::count(celltypes_1, sample_id) %>% # only computing number of cells, not percentages
  reshape2::dcast(sample_id~celltypes_1) %>% # change rows to columns
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) # replacing NA by 0
head(celllineage_prop_df)

write_csv(celllineage_prop_df, file = paste0(path_r_tables_out, "/post12_cd45neg_lineages_counts.csv"))

## Cell-types & sample (treatment)
celltypesprop_df <- seurat_obj_metadata %>% 
  select(c("sample_id", "celltypes_2")) %>%
  dplyr::count(celltypes_2, sample_id) %>% # only computing number of cells, not percentages
  reshape2::dcast(sample_id~celltypes_2) %>% # change rows to columns
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) # replacing NA by 0
head(celltypesprop_df)

write_csv(celltypesprop_df, file = paste0(path_r_tables_out, "/post12_cd45neg_celltypes_counts.csv"))
