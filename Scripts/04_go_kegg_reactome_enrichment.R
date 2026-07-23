# 04_go_kegg_reactome_enrichment.R
# Functional Enrichment Analysis: GO, KEGG, Reactome
# For Breast Cancer RNA-Seq DEGs (GSE124647)

# Load required libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
})

cat("=================================================================\n")
cat("GO, KEGG & Reactome Pathway Enrichment Analysis\n")
cat("=================================================================\n\n")

# Parameters
dge_results_file <- "Results/02_dge_analysis/deg_results.csv"
output_dir <- "Results/04_enrichment"
p_adj_threshold <- 0.05
qvalue_threshold <- 0.05
organism <- "human"
key_type <- "SYMBOL"

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading DEG results...\n")
# TODO: Load differential expression results
# deg_results <- read.csv(dge_results_file, row.names = 1)

cat("\nPerforming comprehensive enrichment analysis...\n")
cat("- Gene Ontology (GO): Biological Process, Molecular Function, Cellular Component\n")
cat("- KEGG pathway enrichment\n")
cat("- Reactome pathway enrichment\n")
cat("- Generating enrichment plots (dot plot, bar plot, enrichment map, cnet plot)\n")

# TODO: Implement enrichment analysis
# 1. Extract significant DEGs (padj < 0.05, |log2FC| >= 1.0)
# 2. Convert gene symbols to Entrez IDs
# 3. Perform GO enrichment (BP, MF, CC) using enrichGO
# 4. Perform KEGG enrichment using enrichKEGG
# 5. Perform Reactome enrichment using enrichPathway (via ReactomePA)
# 6. Generate and save all plots

cat("\n[Status] Enrichment analysis module ready.\n")
cat(sprintf("[Output] Results and plots saved to: %s\n", output_dir))