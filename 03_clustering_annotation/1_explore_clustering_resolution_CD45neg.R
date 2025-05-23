# ---
# author: "Laura Jiménez Gracia"
# date: 2021-04-27
# ---
# This R script allows to group cells into different clusters 
# based on the similarity of their gene expression profiles,
# considering different resolution (granularity) levels.


# Pre-processing
## Load packages
library(tidyverse)
library(Seurat)
library(RColorBrewer)

## Paths
path_r_objects_in <- here::here("02_QC/results/R_objects")
path_r_objects_out <- here::here("03_clustering_annotation/results/R_objects")
path_r_figs <- here::here("03_clustering_annotation/results/figs")

# Functions
source(here::here("bin/utils.R"))

# Parameters
## Define resolutions to compute
resolutions_range <- c(0.01, 0.05, 0.1, 0.2, 0.25,
                       0.3, 0.4, 0.5, 0.6, 0.7, 
                       0.75, 0.8, 0.9, 1, 1.1, 
                       1.2, 1.25, 1.3, 1.4, 1.5)
resolution_names <- as.character(purrr::map(resolutions_range, function(res) {
  paste0("RNA_snn_res.", res)}))

color_palette = c("aquamarine", "orange", "lightgreen", "blueviolet", "maroon",
                  "deepskyblue3", "darkgray", "red", "green", "gold3", "magenta", 
                  "pink", "lightgray", "blue", "darkgreen",  "gold", "coral2", 
                  "pink3", "yellow", "palegreen", "peachpuff", "steelblue", "cyan3",
                  "darkblue", "lightblue", "green3", "aquamarine", "orange",
                  "darkgray", "red", "green", "deepskyblue3",  "gold3", "magenta", 
                  "pink", "violet", "blue", "darkgreen",  "gold", "coral2", 
                  "navy", "yellow", "darkblue", "palegreen", "peachpuff", "steelblue",
                  "cyan3", "lightblue", "green3")

# Load Data
## Load Seurat object
seurat_obj <- readRDS(paste0(path_r_objects_in, "/post12_cd45neg_lognorm_harmony_processed.rds"))
DefaultAssay(seurat_obj) <- "RNA"


# Clustering
## Determine the K-nearest neighbor graph
seurat_obj <- FindNeighbors(
  seurat_obj,
  reduction = "harmony",
  dims = 1:30)

## Determine the clusters (Louvain algorithm) for multiple resolutions                                
seurat_obj <- FindClusters(
  seurat_obj,
  resolution = resolutions_range,
  verbose = FALSE)


# Clustering overview
gg_umap_cluster_resolution <- DimPlot(object = seurat_obj,
                                      group.by = resolution_names,
                                      label = TRUE,
                                      label.size = 3,
                                      cols = color_palette,
                                      ncol = 4
                                      ) & NoLegend()

# Save image
ggsave(filename = paste0(path_r_figs, "/post12_cd45neg_clustering_resolutions_umap.png"),
       plot = gg_umap_cluster_resolution,
       width = 25,
       height = 25)


# Save cell-type biomarkers
## Export in .rds object
saveRDS(seurat_obj, 
        file = paste0(path_r_objects_out, "/post12_cd45neg_clustering_resolutions.rds"))
#seurat_obj <- readRDS(paste0(path_r_objects_out, "/post12_cd45neg_clustering_resolutions.rds"))
