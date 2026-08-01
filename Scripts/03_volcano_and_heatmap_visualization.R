#!/usr/bin/env Rscript
# 03_volcano_and_heatmap_visualization.R
# Module 3: Publication-grade DEG Volcano Plot & Expression Heatmap
# Dataset: GSE124647 Breast Cancer (TNBC vs Luminal)
# Inputs : Data/deseq2_full_results.rds, Data/significant_degs.rds,
#          Data/vst_counts.rds, Data/sample_metadata.rds
# Outputs: Figures/03_volcano_plot.{png,pdf}
#          Figures/03_deg_heatmap.{png,pdf}
# Usage  : Rscript Scripts/03_volcano_and_heatmap_visualization.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(pheatmap)
  library(ggrepel)
  library(RColorBrewer)
  library(scales)
})

# ============================================================================
# GLOBAL PARAMETERS
# ============================================================================
DATA_DIR      <- "Data"
FIG_DIR       <- "Figures"
PADJ_CUTOFF   <- 0.05
LOG2FC_CUTOFF <- 1.5
TOP_LABEL_N   <- 15          # top driver DEGs labeled on volcano
N_HEAT_GENES  <- 100         # top significant DEGs in heatmap

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=================================================================\n")
cat("03_volcano_and_heatmap_visualization.R  |  GSE124647\n")
cat("=================================================================\n\n")

# ============================================================================
# STEP 1: LOAD INPUT DATA (Module 2 / Module 1 outputs)
# ============================================================================
cat("[STEP 1] Loading DESeq2 results + expression matrix...\n")

full_results <- readRDS(file.path(DATA_DIR, "deseq2_full_results.rds"))
sig_degs     <- readRDS(file.path(DATA_DIR, "significant_degs.rds"))
vst_counts   <- readRDS(file.path(DATA_DIR, "vst_counts.rds"))
metadata     <- readRDS(file.path(DATA_DIR, "sample_metadata.rds"))

cat(sprintf("[INFO] Full DESeq2 results : %d genes\n", nrow(full_results)))
cat(sprintf("[INFO] Significant DEGs    : %d genes\n", nrow(sig_degs)))
cat(sprintf("[INFO] VST expression      : %d genes x %d samples\n",
            nrow(vst_counts), ncol(vst_counts)))

# ============================================================================
# STEP 2: VOLCANO PLOT DATA PREPARATION
# ============================================================================
cat("\n[STEP 2] Preparing volcano plot data...\n")

volcano_df <- full_results %>%
  mutate(
    Direction = case_when(
      Padj < PADJ_CUTOFF & Log2FoldChange >=  LOG2FC_CUTOFF ~ "Upregulated",
      Padj < PADJ_CUTOFF & Log2FoldChange <= -LOG2FC_CUTOFF ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    NegLog10Padj = -log10(Padj)
  )

cat(sprintf("[INFO] Upregulated    : %d\n",
            sum(volcano_df$Direction == "Upregulated")))
cat(sprintf("[INFO] Downregulated  : %d\n",
            sum(volcano_df$Direction == "Downregulated")))
cat(sprintf("[INFO] Not significant: %d\n",
            sum(volcano_df$Direction == "Not significant")))

# Top driver genes: most significant DEGs (highest -log10 Padj)
top_drivers <- volcano_df %>%
  filter(Direction != "Not significant") %>%
  arrange(desc(NegLog10Padj)) %>%
  slice_head(n = TOP_LABEL_N) %>%
  pull(Gene)

volcano_df$Label <- ifelse(volcano_df$Gene %in% top_drivers,
                           volcano_df$Gene, "")

cat(sprintf("[INFO] Top %d driver DEGs labeled: %s\n",
            length(top_drivers),
            paste(head(top_drivers, 5), collapse = ", ")))

# ============================================================================
# STEP 3: VOLCANO PLOT (ggplot2 + ggrepel)
# ============================================================================
cat("\n[STEP 3] Generating volcano plot...\n")

# Publication color scheme
volcano_colors <- c("Upregulated"     = "#D7263D",   # vivid red
                    "Downregulated"   = "#1B6CA8",   # deep blue
                    "Not significant" = "#B0B0B0")   # grey

volcano_plot <- ggplot(volcano_df,
                       aes(x = Log2FoldChange, y = NegLog10Padj,
                           color = Direction)) +
  # Non-significant points plotted first (underlay)
  geom_point(data = filter(volcano_df, Direction == "Not significant"),
             size = 1.2, alpha = 0.55) +
  # Significant points overlaid
  geom_point(data = filter(volcano_df, Direction != "Not significant"),
             size = 2.0, alpha = 0.85) +
  # Cutoff reference lines
  geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(PADJ_CUTOFF),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  # Driver gene labels
  geom_text_repel(aes(label = Label),
                  size = 3.2, max.overlaps = 20,
                  segment.color = "grey50", segment.size = 0.3,
                  box.padding = 0.4, force = 2,
                  show.legend = FALSE) +
  scale_color_manual(values = volcano_colors) +
  scale_x_continuous(limits = c(-max(abs(volcano_df$Log2FoldChange)) * 1.05,
                                max(abs(volcano_df$Log2FoldChange)) * 1.05)) +
  labs(
    title    = "Volcano Plot: TNBC vs Luminal (GSE124647)",
    subtitle = sprintf("Cutoffs: padj < %.2f  |  |log2FC| >= %.1f  |  Top %d driver DEGs",
                       PADJ_CUTOFF, LOG2FC_CUTOFF, TOP_LABEL_N),
    x        = expression(log[2]~Fold~Change ~ (TNBC / Luminal)),
    y        = expression(-log[10]~adjusted~p-value),
    color    = "Regulation"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "03_volcano_plot.png"),
       plot = volcano_plot, width = 10, height = 8, dpi = 300)
pdf(file.path(FIG_DIR, "03_volcano_plot.pdf"),
    width = 10, height = 8)
print(volcano_plot)
dev.off()

cat("[OUTPUT] Figures/03_volcano_plot.png\n")
cat("[OUTPUT] Figures/03_volcano_plot.pdf\n")

# ============================================================================
# STEP 4: HEATMAP DATA PREPARATION (TOP DEGs, Z-SCORED VST)
# ============================================================================
cat("\n[STEP 4] Preparing heatmap matrix (top DEGs, z-score)...\n")

# Select top significant DEGs ranked by adjusted p-value
top_degs <- sig_degs %>%
  arrange(Padj) %>%
  slice_head(n = N_HEAT_GENES)

heat_genes <- intersect(top_degs$Gene, rownames(vst_counts))
if (length(heat_genes) < 10) {
  stop("Fewer than 10 DEGs matched VST matrix rownames; check identifier alignment")
}
cat(sprintf("[INFO] Matched %d DEGs in VST matrix\n", length(heat_genes)))

heat_matrix <- vst_counts[heat_genes, , drop = FALSE]

# Row-wise z-score standardization for maximum contrast
heat_mat_z <- t(scale(t(heat_matrix)))
heat_mat_z <- heat_mat_z[complete.cases(heat_mat_z), , drop = FALSE]

# ============================================================================
# STEP 5: CLINICAL ANNOTATION (TNBC vs Luminal)
# ============================================================================
cat("\n[STEP 5] Adding clinical subtype annotations...\n")

# Align sample annotation to heatmap columns
sample_names <- colnames(heat_mat_z)
if ("source_name" %in% colnames(metadata)) {
  cond_vec <- metadata$source_name[match(sample_names, rownames(metadata))]
} else {
  cond_vec <- metadata$condition[match(sample_names, rownames(metadata))]
}
cond_vec[is.na(cond_vec)] <- "Unknown"

ann_col <- data.frame(
  Subtype = factor(cond_vec),
  row.names = sample_names
)

ann_colors <- list(
  Subtype = c("TNBC" = "#E64B35",
              "Luminal" = "#00A087",
              "Unknown" = "#BBBBBB")
)
cat(sprintf("[INFO] Heatmap columns: %d samples annotated\n", nrow(ann_col)))

# ============================================================================
# STEP 6: HEATMAP (pheatmap, hierarchical clustering)
# ============================================================================
cat("\n[STEP 6] Generating hierarchical clustering heatmap...\n")

heat_colors <- colorRampPalette(c("#3B6FB5", "white", "#FF6B6B"))(100)

# PNG
pheatmap(heat_mat_z,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = heat_colors,
         breaks = seq(-3, 3, length.out = 101),
         fontsize_row = 6,
         fontsize_col = 1,
         border_color = NA,
         main = "Top 100 DEG Expression Heatmap (Z-score VST)",
         annotation_legend = TRUE,
         legend = TRUE,
         width = 10, height = 12,
         filename = file.path(FIG_DIR, "03_deg_heatmap.png"))

# PDF
pheatmap(heat_mat_z,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = heat_colors,
         breaks = seq(-3, 3, length.out = 101),
         fontsize_row = 6,
         fontsize_col = 1,
         border_color = NA,
         main = "Top 100 DEG Expression Heatmap (Z-score VST)",
         annotation_legend = TRUE,
         legend = TRUE,
         width = 10, height = 12,
         filename = file.path(FIG_DIR, "03_deg_heatmap.pdf"))

cat("[OUTPUT] Figures/03_deg_heatmap.png\n")
cat("[OUTPUT] Figures/03_deg_heatmap.pdf\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n=================================================================\n")
cat("VISUALIZATION SUMMARY\n")
cat("=================================================================\n")
cat(sprintf("Volcano plot      : %d up / %d down / %d NS\n",
            sum(volcano_df$Direction == "Upregulated"),
            sum(volcano_df$Direction == "Downregulated"),
            sum(volcano_df$Direction == "Not significant")))
cat(sprintf("Top drivers labeled: %d\n", length(top_drivers)))
cat(sprintf("Heatmap genes     : %d significant DEGs (z-scored)\n", nrow(heat_mat_z)))
cat(sprintf("Heatmap annotations: %s (%d) vs %s (%d)\n",
            levels(ann_col$Subtype)[1],
            sum(ann_col$Subtype == levels(ann_col$Subtype)[1]),
            levels(ann_col$Subtype)[2],
            sum(ann_col$Subtype == levels(ann_col$Subtype)[2])))
cat("\n[SUCCESS] Module 3 visualizations generated!\n")
cat("=================================================================\n")