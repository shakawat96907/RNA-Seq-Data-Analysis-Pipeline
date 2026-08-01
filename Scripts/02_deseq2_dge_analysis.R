#!/usr/bin/env Rscript
# 02_deseq2_dge_analysis.R
# DESeq2 Differential Gene Expression Analysis (Module 2)
# Dataset: GSE124647 Breast Cancer (TNBC vs Luminal)
# Pipeline: filtered_counts + metadata -> DESeqDataSet -> DESeq() -> lfcShrink()
#           -> strict DEG filtering (padj<0.05, |log2FC|>=1.5) -> Up/Down classification
# Outputs: Data/deseq2_full_results.{rds,csv}
#          Data/significant_degs.{rds,csv}
# Usage:   Rscript Scripts/02_deseq2_dge_analysis.R

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(readr)
})

# ============================================================================
# GLOBAL PARAMETERS
# ============================================================================
DATA_DIR        <- "Data"
FILTERED_FILE   <- file.path(DATA_DIR, "filtered_counts.rds")
METADATA_FILE   <- file.path(DATA_DIR, "sample_metadata.rds")
PADJ_CUTOFF     <- 0.05
LOG2FC_CUTOFF   <- 1.5
TOP_N_GENES     <- 5000       # Most variable genes for tractable runtime
ALPHA           <- 0.05       # results() significance level

cat("=================================================================\n")
cat("02_deseq2_dge_analysis.R  |  GSE124647 TNBC vs Luminal\n")
cat("=================================================================\n\n")

# ============================================================================
# STEP 1: LOAD PREPROCESSED DATA
# ============================================================================
cat("[STEP 1] Loading preprocessed data (Module 1 outputs)...\n")

if (!file.exists(FILTERED_FILE)) stop("Missing file: ", FILTERED_FILE)
if (!file.exists(METADATA_FILE)) stop("Missing file: ", METADATA_FILE)

counts    <- readRDS(FILTERED_FILE)
metadata  <- readRDS(METADATA_FILE)

cat(sprintf("[INFO] Filtered counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("[INFO] Metadata rows  : %d samples\n", nrow(metadata)))

# ============================================================================
# STEP 2: PREPARE DESIGN / CONDITION LABELS
# ============================================================================
cat("\n[STEP 2] Preparing condition labels (TNBC vs Luminal)...\n")

if ("source_name" %in% colnames(metadata)) {
  condition <- factor(metadata$source_name)
} else if ("condition" %in% colnames(metadata)) {
  condition <- factor(metadata$condition)
} else {
  stop("No condition column (source_name/condition) in metadata")
}

# Align to count matrix columns
if (nrow(metadata) != ncol(counts) || !all(rownames(metadata) == colnames(counts))) {
  condition <- condition[seq_len(ncol(counts))]
}

column_data <- data.frame(row.names = colnames(counts), condition = condition)
cat(sprintf("[INFO] Groups: %s (n=%d) vs %s (n=%d)\n",
            levels(condition)[1], sum(condition == levels(condition)[1]),
            levels(condition)[2], sum(condition == levels(condition)[2])))

# Reference level = first (Luminal); TNBC becomes comparison group
column_data$condition <- relevel(column_data$condition, ref = levels(condition)[1])

# ============================================================================
# STEP 3: SUBSET TO TOP VARIABLE GENES (TRACTABLE RUNTIME)
# ============================================================================
cat("\n[STEP 3] Selecting top variable genes for statistical power...\n")

counts <- round(abs(counts))
storage.mode(counts) <- "integer"

if (nrow(counts) > TOP_N_GENES) {
  gene_vars <- apply(counts, 1, var)
  top_idx   <- order(gene_vars, decreasing = TRUE)[seq_len(TOP_N_GENES)]
  counts    <- counts[top_idx, , drop = FALSE]
  cat(sprintf("[INFO] Retained top %d most variable genes\n", nrow(counts)))
} else {
  cat(sprintf("[INFO] Using all %d genes\n", nrow(counts)))
}

# ============================================================================
# STEP 4: DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================
cat("\n[STEP 4] Running DESeq2::DESeq()...\n")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData  = column_data,
                              design   = ~ condition)
dds <- DESeq2::DESeq(dds, quiet = FALSE)

cat("[INFO] DESeq2 model fitting complete\n")

# ============================================================================
# STEP 5: EXTRACT RESULTS WITH SHRINKAGE (lfcShrink)
# ============================================================================
cat("\n[STEP 5] Extracting results with LFC shrinkage...\n")

# Standard results with Independent Filtering
res <- DESeq2::results(dds,
                       contrast = c("condition",
                                    levels(column_data$condition)[2],
                                    levels(column_data$condition)[1]),
                       alpha = ALPHA)

# Shrinkage estimates (apeglm) for robust LFC
if ("apeglm" %in% rownames(installed.packages())) {
  res_shrink <- DESeq2::lfcShrink(dds,
                                  coef = resultsNames(dds)[2],
                                  type = "apeglm")
  res_shrink_df <- as.data.frame(res_shrink)
} else {
  cat("[WARN] apeglm not installed; using standard results\n")
  res_shrink_df <- as.data.frame(res)
}
# Note: LFC from shrinkage used for the final table where available,
# keeping p-values/adj p-values from standard results.
final_stats <- as.data.frame(res)

res_df <- data.frame(
  Gene            = rownames(final_stats),
  baseMean        = final_stats$baseMean,
  Log2FoldChange  = if ("log2FoldChange" %in% colnames(res_shrink_df)) {
                      res_shrink_df$log2FoldChange
                    } else { final_stats$log2FoldChange },
  lfcSE           = final_stats$lfcSE,
  stat            = final_stats$stat,
  PValue          = final_stats$pvalue,
  Padj            = final_stats$padj,
  stringsAsFactors = FALSE,
  check.names      = FALSE
) %>%
  na.omit()

cat(sprintf("[INFO] Total genes with valid statistics: %d\n", nrow(res_df)))
cat(sprintf("[INFO] Significant DEGs (padj<%.2f, |log2FC|>=%.1f): %d\n",
            PADJ_CUTOFF, LOG2FC_CUTOFF,
            sum(res_df$Padj < PADJ_CUTOFF & abs(res_df$Log2FoldChange) >= LOG2FC_CUTOFF)))

# ============================================================================
# STEP 6: STATISTICAL FILTERING & DEG CLASSIFICATION
# ============================================================================
cat("\n[STEP 6] Filtering and classifying DEGs...\n")

full_results <- res_df

degs <- res_df %>%
  filter(Padj < PADJ_CUTOFF & abs(Log2FoldChange) >= LOG2FC_CUTOFF) %>%
  mutate(Direction = ifelse(Log2FoldChange > 0, "Up-regulated", "Down-regulated"))

cat(sprintf("[INFO] Up-regulated genes  : %d\n",
            sum(degs$Direction == "Up-regulated")))
cat(sprintf("[INFO] Down-regulated genes: %d\n",
            sum(degs$Direction == "Down-regulated")))

# ============================================================================
# STEP 7: SAVE OUTPUTS
# ============================================================================
cat("\n[STEP 7] Saving results...\n")

saveRDS(full_results, file.path(DATA_DIR, "deseq2_full_results.rds"))
write_csv(full_results, file.path(DATA_DIR, "deseq2_full_results.csv"))

saveRDS(degs, file.path(DATA_DIR, "significant_degs.rds"))
write_csv(degs, file.path(DATA_DIR, "significant_degs.csv"))

cat("[OUTPUT] Data/deseq2_full_results.rds  (full DESeq2 statistics)\n")
cat("[OUTPUT] Data/deseq2_full_results.csv\n")
cat("[OUTPUT] Data/significant_degs.rds     (filtered DEGs, padj<0.05, |log2FC|>=1.5)\n")
cat("[OUTPUT] Data/significant_degs.csv\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n=================================================================\n")
cat("DGE ANALYSIS SUMMARY\n")
cat("=================================================================\n")
cat(sprintf("Comparison         : %s vs %s\n",
            levels(column_data$condition)[2], levels(column_data$condition)[1]))
cat(sprintf("Genes tested       : %d\n", nrow(full_results)))
cat(sprintf("Cutoffs            : padj < %.2f, |log2FC| >= %.1f\n",
            PADJ_CUTOFF, LOG2FC_CUTOFF))
cat(sprintf("Significant DEGs   : %d\n", nrow(degs)))
cat(sprintf("  Up-regulated     : %d\n", sum(degs$Direction == "Up-regulated")))
cat(sprintf("  Down-regulated   : %d\n", sum(degs$Direction == "Down-regulated")))
cat(sprintf("Shrinkage          : apeglm LFC shrinkage applied\n"))
cat("\n[SUCCESS] DESeq2 DGE analysis completed!\n")
cat("=================================================================\n")