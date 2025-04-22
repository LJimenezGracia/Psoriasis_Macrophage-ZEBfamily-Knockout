# Psoriasis_Macrophage-ZEBfamily-Knockout

This is the GitHub repository for the study to understand the **Impact of Zeb1/2 Knock-out in Macrophages on Psoriasis Pathogenesis**.

The manuscript is currently under a review process.

## Abstract


## Code implementation

The repository is organized into the following folder tree, which contains all the necessary files and scripts to perform the detailed tasks and reproduce all our results

* **01_cellranger_mapping** --> It includes an overview of the project data information, including which samples and 10X libraries were generated. Following the 10X Genomics CellPlex strategy, here we also obtain the demultiplex count matrix for each sample and library. Also, it contains the scripts needed to create a folder directory to perform the sequencing read mapping to the reference genome. 

* **02_QC** --> R markdown notebooks to perform the quality control and the first pre-processing, including data normalization, scaling, dimensionality reduction and integration.

* **03_clustering_annotation** --> All R markdown notebooks to perform a top-down clustering annotation approach, as well as scripts to find differential expressed markers for each clustering, and to assign a biological-relevant identity to each cluster.

* **04_analysis** --> All the code used to perform further downstream analysis on the GEX processed data. It includes all the scripts and R markdown notebooks to perform cell composition analysis, differential expression analysis (DEA) followed by gene set enrichment analysis (GSEA), as well as gene expression and gene signature evaluation. 


### Package versions

The (most important) packages and versions needed to reproduce the full analysis are listed below:

* CellRanger (v6.1.1) was used to mapped single-cell RNA-seq reads (10X Genomics) to the reference genome.

*--- in R (v4.0.5) ---*
* Seurat (> v4.0.0)
* SeuratObject (v4.0.1)
* Harmony (v1.0)
* fgsea (v1.16)
* Ucell (v2.2.0)
* ggplot2 (v3.3.3)

*--- in Python (v3.9.19) ---*
* scanpy (v1.10.3)
* scCODA (v0.1.9)

## Data accessibility

* The complete raw data (FASTQ files) generated in this study, as well as the processed count matrices, have been submitted to the NCBI Gene Expression Omnibus (GEO) under accession number [XXXX](XXXX).


## Code accessibility

You can easily download a copy of all the files contained in this repository by:

* Cloning the git repository using the following command in the terminal:

`git clone https://github.com/LJimenezGracia/Psoriasis_Macrophage-ZEBfamily-Knockout.git`

* Downloading a .ZIP archive [HERE](https://github.com/LJimenezGracia/Psoriasis_Macrophage-ZEBfamily-Knockout/archive/refs/heads/main.zip).
