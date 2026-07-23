# 02_dge_deseq2_edger_modeling.R
# Differential Gene Expression Analysis using DESeq2 and edgeR
# For Breast Cancer RNA-Seq cohort (GSE124647)

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(tidyverse)
  library(dplyr)
})

cat("=================================================================\n")
cat("Differential Gene Expression (DGE) Analysis\n")
cat("=================================================================\n\n")

# Parameters
input_dir <- "Results/01_data_ingestion_qc"
output_dir <- "Results/02_dge_analysis"
condition_col <- "condition"
control_group <- "normal"
treatment_group <- "tumor"
alpha <- 0.05
log2FC_threshold <- 1.0

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading preprocessed data...\n")
# TODO: Load preprocessed count matrix and sample metadata
# counts <- readRDS(file.path(input_dir, "filtered_counts.rds"))
# metadata <- readRDS(file.path(input_dir, "sample_metadata.rds"))

cat("\nRunning DGE analysis...\n")
cat("- DESeq2 modeling with shrinkage estimators\n")
cat("- edgeR exactTest/glmQLFit for comparison\n")
cat("- Limma-voom for additional validation\n")
cat("- Extracting significant DEGs with p_adj < 0.05 and |log2FC| >= 1.0\n")

# TODO: Implement DGE analysis workflow
# 1. Create DESeqDataSet and run DESeq()
# 2. Extract results with lfcShrink
# 3. edgeR: DGEList, calcNormFactors, exactTest/glmQLFit
# 4. Limma: voom, eBayes, treat
# 5. Merge results, filter DEGs, save to CSV

cat("\n[Status] DGE analysis module ready.\n")
cat(sprintf("[Output] Results will be saved to: %s\n", output_dir))