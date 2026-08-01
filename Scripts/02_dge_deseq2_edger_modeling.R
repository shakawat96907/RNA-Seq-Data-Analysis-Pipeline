#!/usr/bin/env Rscript
# 02_dge_deseq2_edger_modeling.R
# Differential Gene Expression Analysis using DESeq2 and edgeR
# Breast Cancer GSE124647: TNBC vs Luminal

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

cat("=================================================================\n")
cat("DGE Analysis: DESeq2 + edgeR Consensus\n")
cat("=================================================================\n\n")

# ============================================================================
# PARAMETERS
# ============================================================================
raw_counts_file <- "Data/raw_counts.rds"
metadata_file <- "Data/sample_metadata.rds"
output_dir <- "Data"

padj_threshold <- 0.05
log2fc_threshold <- 1.5

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================
cat("[STEP 1] Loading preprocessed data...\n")

if (!file.exists(raw_counts_file)) {
  stop("File not found: ", raw_counts_file)
}
if (!file.exists(metadata_file)) {
  stop("File not found: ", metadata_file)
}

counts <- readRDS(raw_counts_file)
metadata <- readRDS(metadata_file)

cat(sprintf("[INFO] Count matrix: %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("[INFO] Metadata: %d samples\n", nrow(metadata)))

# ============================================================================
# STEP 2: ENSURE INTEGER COUNTS AND PREPARE DESIGN MATRIX
# ============================================================================
cat("\n[STEP 2] Preparing design matrix...\n")

# Ensure counts are non-negative integers (required by DESeq2 and edgeR)
counts <- round(abs(counts))
storage.mode(counts) <- "integer"
cat(sprintf("[INFO] Counts converted to integer type (min: %d, max: %d)\n", min(counts), max(counts)))

# Extract condition labels
if ("source_name" %in% colnames(metadata)) {
  condition <- metadata$source_name
} else if ("condition" %in% colnames(metadata)) {
  condition <- metadata$condition
} else {
  stop("No condition column found in metadata")
}

condition <- factor(condition)
cat(sprintf("[INFO] Conditions: %s\n", paste(levels(condition), collapse = " vs ")))

# Build design matrix
design <- model.matrix(~ condition)
colnames(design) <- gsub("condition", "", colnames(design))

# ============================================================================
# STEP 3: SUBSET TO TOP VARIABLE GENES FOR TRACTABLE RUNTIME
# ============================================================================
cat("\n[STEP 3] Subsetting to top variable genes for tractable runtime...\n")

# Select top 5000 most variable genes to keep runtime manageable
if (nrow(counts) > 5000) {
  gene_vars <- apply(counts, 1, var)
  top_genes <- order(gene_vars, decreasing = TRUE)[1:5000]
  counts <- counts[top_genes, , drop = FALSE]
  cat(sprintf("[INFO] Subset to %d most variable genes\n", nrow(counts)))
} else {
  cat(sprintf("[INFO] Using all %d genes\n", nrow(counts)))
}

# ============================================================================
# STEP 4: DESEQ2 ANALYSIS
# ============================================================================
cat("\n[STEP 4] Running DESeq2...\n")

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = data.frame(condition = condition),
                              design = ~ condition)

# Run DESeq2 pipeline
dds <- DESeq2::DESeq(dds)

# Extract results
res <- DESeq2::results(dds, 
                       contrast = c("condition", levels(condition)[2], levels(condition)[1]),
                       independentFiltering = TRUE)

# Standardize column names
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene") %>%
  mutate(Log2FoldChange = log2FoldChange,
         PValue = pvalue,
         Padj = padj) %>%
  select(Gene, Log2FoldChange, PValue, Padj)

cat(sprintf("[DESeq2] Total genes tested: %d\n", nrow(res_df)))
cat(sprintf("[DESeq2] DEGs (padj < %.2f, |log2FC| >= %.1f): %d\n", 
            padj_threshold, log2fc_threshold, 
            sum(res_df$Padj < padj_threshold & abs(res_df$Log2FoldChange) >= log2fc_threshold)))

# Save DESeq2 results
saveRDS(res_df, file.path(output_dir, "deseq2_results.rds"))
write_csv(res_df, file.path(output_dir, "deseq2_results.csv"))
cat("[OUTPUT] Saved: Data/deseq2_results.rds, Data/deseq2_results.csv\n")

# ============================================================================
# STEP 5: EDGER ANALYSIS
# ============================================================================
cat("\n[STEP 5] Running edgeR...\n")

# Create DGEList
dge <- DGEList(counts = counts, group = condition)
dge <- calcNormFactors(dge, method = "TMM")

# Estimate dispersion
dge <- estimateDisp(dge, design)

# glmQLFit
fit <- glmQLFit(dge, design)
ql_test <- glmQLFTest(fit, coef = 2)

# Extract results
edger_df <- topTags(ql_test, n = Inf)$table %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  mutate(Log2FoldChange = logFC,
         PValue = PValue,
         Padj = FDR) %>%
  select(Gene, Log2FoldChange, PValue, Padj)

cat(sprintf("[edgeR] Total genes tested: %d\n", nrow(edger_df)))
cat(sprintf("[edgeR] DEGs (padj < %.2f, |log2FC| >= %.1f): %d\n", 
            padj_threshold, log2fc_threshold,
            sum(edger_df$Padj < padj_threshold & abs(edger_df$Log2FoldChange) >= log2fc_threshold)))

# Save edgeR results
saveRDS(edger_df, file.path(output_dir, "edger_results.rds"))
write_csv(edger_df, file.path(output_dir, "edger_results.csv"))
cat("[OUTPUT] Saved: Data/edger_results.rds, Data/edger_results.csv\n")

# ============================================================================
# STEP 6: CONSENSUS DEG IDENTIFICATION
# ============================================================================
cat("\n[STEP 6] Identifying consensus DEGs...\n")

# Filter significant DEGs for each method
deseq_sig <- res_df %>%
  filter(Padj < padj_threshold & abs(Log2FoldChange) >= log2fc_threshold) %>%
  mutate(Direction = ifelse(Log2FoldChange > 0, "Up-regulated", "Down-regulated"))

edger_sig <- edger_df %>%
  filter(Padj < padj_threshold & abs(Log2FoldChange) >= log2fc_threshold) %>%
  mutate(Direction = ifelse(Log2FoldChange > 0, "Up-regulated", "Down-regulated"))

cat(sprintf("[DESeq2] Significant DEGs: %d (Up: %d, Down: %d)\n",
            nrow(deseq_sig),
            sum(deseq_sig$Direction == "Up-regulated"),
            sum(deseq_sig$Direction == "Down-regulated")))

cat(sprintf("[edgeR] Significant DEGs: %d (Up: %d, Down: %d)\n",
            nrow(edger_sig),
            sum(edger_sig$Direction == "Up-regulated"),
            sum(edger_sig$Direction == "Down-regulated")))

# Consensus: genes significant in BOTH methods with same directionality
consensus_up <- intersect(deseq_sig$Gene[deseq_sig$Direction == "Up-regulated"],
                          edger_sig$Gene[edger_sig$Direction == "Up-regulated"])
consensus_down <- intersect(deseq_sig$Gene[deseq_sig$Direction == "Down-regulated"],
                            edger_sig$Gene[edger_sig$Direction == "Down-regulated"])
consensus_genes <- union(consensus_up, consensus_down)

# Build consensus table from DESeq2 results
consensus <- res_df %>%
  filter(Gene %in% consensus_genes) %>%
  mutate(Direction = ifelse(Log2FoldChange > 0, "Up-regulated", "Down-regulated"))

cat(sprintf("[Consensus] DEGs in both methods: %d\n", nrow(consensus)))
cat(sprintf("[Consensus] Up-regulated: %d\n", sum(consensus$Direction == "Up-regulated")))
cat(sprintf("[Consensus] Down-regulated: %d\n", sum(consensus$Direction == "Down-regulated")))

# Save consensus results
saveRDS(consensus, file.path(output_dir, "consensus_degs.rds"))
write_csv(consensus, file.path(output_dir, "consensus_degs.csv"))
cat("[OUTPUT] Saved: Data/consensus_degs.rds, Data/consensus_degs.csv\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n=================================================================\n")
cat("DGE ANALYSIS SUMMARY\n")
cat("=================================================================\n")
cat(sprintf("Comparison:               %s vs %s\n", levels(condition)[2], levels(condition)[1]))
cat(sprintf("Total genes tested:       %d\n", nrow(res_df)))
cat(sprintf("DESeq2 DEGs:              %d\n", nrow(deseq_sig)))
cat(sprintf("edgeR DEGs:               %d\n", nrow(edger_sig)))
cat(sprintf("Consensus DEGs:            %d\n", nrow(consensus)))
cat(sprintf("  Up-regulated:           %d\n", sum(consensus$Direction == "Up-regulated")))
cat(sprintf("  Down-regulated:         %d\n", sum(consensus$Direction == "Down-regulated")))
cat(sprintf("Cutoffs:                  padj < %.2f, |log2FC| >= %.1f\n", padj_threshold, log2fc_threshold))
cat("\n[SUCCESS] DGE analysis completed!\n")
cat("=================================================================\n")