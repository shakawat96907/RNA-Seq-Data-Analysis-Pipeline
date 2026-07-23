# 01_data_ingestion_qc_preprocessing.R
# Data Ingestion, Quality Control, and Preprocessing for Breast Cancer RNA-Seq
# Target cohort: GSE124647

# Load required libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(tidyverse)
  library(pheatmap)
  library(gridExtra)
})

cat("=================================================================\n")
cat("Data Ingestion, QC & Preprocessing\n")
cat("=================================================================\n\n")

# Parameters
geo_accession <- "GSE124647"
output_dir <- "Results/01_data_ingestion_qc"
counts_dir <- "Data/counts"
metadata_dir <- "Data/metadata"
n_top_genes <- 500

# Create output directories
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Target GEO dataset: %s\n", geo_accession))
cat("\nPipeline stages:\n")
cat("1. Download GEO series matrix and supplementary files\n")
cat("2. Extract count matrix and sample annotations\n")
cat("3. Perform quality control: filtering, normalization, outlier detection\n")
cat("4. Generate sample-level QC metrics and visualizations\n")

# TODO: Implement data ingestion and QC
# 1. getGEO(geo_accession) to fetch metadata
# 2. download supplementary count files
# 3. Create DGEList, filter low-count genes
# 4. Calculate library sizes, TMM normalization
# 5. Generate PCA, MDS, and sample clustering QC plots

cat("\n[Status] Data ingestion and preprocessing module ready.\n")
cat(sprintf("[Output] Results will be saved to: %s\n", output_dir))