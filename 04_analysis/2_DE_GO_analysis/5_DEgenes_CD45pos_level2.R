# ---
# author: "Laura Jiménez Gracia"
# date: 2021-12-04
# ---
# This R script allows to perform a Differential Expression Analysis (DEA)
# between the different groups associated to a particular condition, and
# to find transcriptional signatures associated to that condition.


# Pre-processing
## Load packages
library(tidyverse)
library(Seurat)
library(ggrepel)

## Parameters
# Paths
path_r_objects_in <- here::here("03_clustering_annotation/results/R_objects")
path_r_tables <- here::here("04_analysis/2_DE_GO_analysis/results/tables")
path_r_objects_out <- here::here("04_analysis/2_DE_GO_analysis/results/R_objects")

# Functions
source(here::here("bin/utils.R"))

# Define comparison parameters
comparison_name_list <- c("validate-genotype", "validate-treatment", "compare-IMQKO-Z1Z2",
			"compare-Z1-treatment", "compare-Z1-KO", "compare-Z2-treatment", "compare-Z2-KO")
group_1_list <- c("Z1_WT_CD45pos", "Z1_WT_IMQ_CD45pos", "Z1_KO_IMQ_CD45pos",
			"Z1_WT_IMQ_CD45pos", "Z1_KO_IMQ_CD45pos", "Z2_WT_IMQ_CD45pos", "Z2_KO_IMQ_CD45pos")
group_2_list <- c("Z2_WT_CD45pos", "Z2_WT_IMQ_CD45pos", "Z2_KO_IMQ_CD45pos",
			"Z1_WT_CD45pos", "Z1_WT_IMQ_CD45pos", "Z2_WT_CD45pos", "Z2_WT_IMQ_CD45pos")

## Load data -- not computing for doublets
seurat_obj <- readRDS(paste0(path_r_objects_in, "/post12_cd45pos_alldata_cleaned_final_clustering_annotation.rds"))


for (index in 1:length(comparison_name_list)) {
  comparison_name <- comparison_name_list[[index]]
  group_1 <- group_1_list[[index]]
  group_2 <- group_2_list[[index]]
  
  # Define comparison paths
  path_object <- paste0(path_r_objects_out, "/post12_cd45pos_DEgenes_celltypes_", comparison_name, ".rds")
  path_file <- paste0(path_r_tables, "/post12_cd45pos_DEgenes_celltypes_", comparison_name, ".xlsx")
  path_file_split <- paste0(path_r_tables, "/cd45pos_celltypes/post12_cd45pos_DEgenes_celltypes_", comparison_name, "_UPDOWN.xlsx")
  
  # DE analysis considering cell-types separately
  ## Subset data
  cells_to_keep <- (seurat_obj$sample_id == group_1 | seurat_obj$sample_id == group_2)
  seurat_obj_subset <- seurat_obj[, cells_to_keep]
  
  DE_genes_bycell <- purrr::map(levels(seurat_obj_subset$celltypes_2), function(celltype) {
    seurat_obj_subset <- seurat_obj_subset[, seurat_obj_subset$celltypes_2 == celltype]    
    Idents(seurat_obj_subset) <- "sample_id"
    
    ## Find DE genes
    DE_genes <- find_DEgenes_celltypes(
      seurat = seurat_obj_subset,
      ident_1 = group_1,
      ident_2 = group_2,
      test_de = "wilcox",
      threshold_pvaladj = 0.05,
      min_cells = 0 # default 3, reducing as some conditions up to 1
    )
    
    })
  names(DE_genes_bycell) <- levels(seurat_obj_subset$celltypes_2)
  
  # Save results
  ## DE genes as RDS object
  saveRDS(DE_genes_bycell, file = path_object)
  
  ## Export list of DE genes
  openxlsx::write.xlsx(DE_genes_bycell, file = path_file)


  
  # ## SPLIT up- & down-regulated
  # DE_genes_split <- purrr::map(names(DE_genes_bycell), function(celltype) {
  #   # DEA
  #   DE_genes_up <- DE_genes_bycell[[celltype]] %>% 
  #     dplyr::filter(is_significant & avg_log2FC > 0.25) %>% 
  #     dplyr::select(gene, avg_log2FC, p_val_adj)
  #   
  #   DE_genes_down <- DE_genes_bycell[[celltype]] %>% 
  #     dplyr::filter(is_significant & avg_log2FC < -0.25)%>% 
  #     dplyr::select(gene, avg_log2FC, p_val_adj)
  #   
  #   DE_genes_list <- list(DE_genes_up, DE_genes_down)
  #   names(DE_genes_list) <- c("UP", "DOWN")
  #   return(DE_genes_list)
  #   }
  #   )
  # names(DE_genes_split) <- names(DE_genes_bycell)
  # DE_genes_split_list <- do.call(c, DE_genes_split)
  # DE_genes_split_list <- DE_genes_split_list %>%
  #   discard(is.null)
  # 
  # # Save results
  # ## Export list of DE genes
  # openxlsx::write.xlsx(DE_genes_split_list, file = path_file_split)
}


