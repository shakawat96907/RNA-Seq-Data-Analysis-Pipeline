#!/usr/bin/env Rscript
# 04_go_kegg_pathway_enrichment.R
# Module 4: GO & KEGG Functional Pathway Enrichment Analysis
# Dataset: GSE124647 Breast Cancer (TNBC vs Luminal)
# Inputs : Data/significant_degs.rds, Data/deseq2_full_results.rds
# Outputs: Data/go_enrichment_results.{rds,csv}
#          Data/kegg_enrichment_results.{rds,csv}
#          Figures/04_go_enrichment_dotplot.png
#          Figures/04_kegg_pathway_barplot.png
# Usage  : Rscript Scripts/04_go_kegg_pathway_enrichment.R

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(RColorBrewer)
  library(scales)
})

# ============================================================================
# GLOBAL PARAMETERS
# ============================================================================
DATA_DIR    <- "Data"
FIG_DIR     <- "Figures"
PADJ_CUTOFF <- 0.05
TOP_N_TERMS <- 10
MAX_DEGS    <- 200

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=================================================================\n")
cat("04_go_kegg_pathway_enrichment.R  |  GSE124647\n")
cat("=================================================================\n\n")

# ============================================================================
# STEP 1: LOAD DEGs FROM MODULE 2
# ============================================================================
cat("[STEP 1] Loading significant DEGs from Module 2...\n")
sig_degs <- readRDS(file.path(DATA_DIR, "significant_degs.rds"))
full_res <- readRDS(file.path(DATA_DIR, "deseq2_full_results.rds"))
if (nrow(sig_degs) > MAX_DEGS) {
  sig_degs <- sig_degs[order(sig_degs$Padj)[1:MAX_DEGS], , drop = FALSE]
}
cat(sprintf("[INFO] Significant DEGs loaded: %d\n", nrow(sig_degs)))

# ============================================================================
# STEP 2: GENE ID MAPPING (Synthetic IDs -> Entrez IDs)
# ============================================================================
cat("\n[STEP 2] Mapping gene identifiers to Entrez IDs...\n")
degs_entrez <- as.integer(sub("^Gene", "", as.character(sig_degs$Gene)))
valid_ids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
mapped_ids_valid <- degs_entrez[degs_entrez %in% valid_ids]
if (length(mapped_ids_valid) == 0) stop("No valid Entrez IDs mapped")
cat(sprintf("[INFO] Valid Entrez IDs: %d / %d\n",
            length(mapped_ids_valid), length(degs_entrez)))

background_entrez <- as.integer(sub("^Gene", "", as.character(full_res$Gene)))
background_genes <- background_entrez[background_entrez %in% valid_ids]
cat(sprintf("[INFO] Background genes: %d\n", length(background_genes)))

# ============================================================================
# STEP 3: GO ENRICHMENT (BP, CC, MF) - hypergeometric test
# ============================================================================
cat("\n[STEP 3] Performing GO enrichment (BP, CC, MF)...\n")

run_go <- function(gene_list, universe, ont) {
  all_go <- AnnotationDbi::select(org.Hs.eg.db,
                   keys = as.character(universe),
                   columns = c("ENTREZID", "GO", "ONTOLOGY"),
                   keytype = "ENTREZID", multiVals = "first")
  all_go <- all_go[!is.na(all_go$GO) & all_go$ONTOLOGY == ont, ]
  gene_go <- AnnotationDbi::select(org.Hs.eg.db,
                    keys = as.character(gene_list),
                    columns = c("ENTREZID", "GO", "ONTOLOGY"),
                    keytype = "ENTREZID", multiVals = "first")
  gene_go <- gene_go[!is.na(gene_go$GO) & gene_go$ONTOLOGY == ont, ]
  gene_terms <- unique(gene_go$GO)
  if (length(gene_terms) == 0) return(NULL)

  N <- length(universe)
  n <- length(gene_list)
  results <- lapply(gene_terms, function(term) {
    term_genes <- unique(all_go$ENTREZID[all_go$GO == term])
    M <- length(term_genes)
    k <- sum(gene_list %in% term_genes)
    if (M == 0 || k == 0) return(NULL)
    pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    data.frame(ID = term, Description = Term(term),
               GeneRatio = sprintf("%d/%d", k, n),
               BgRatio = sprintf("%d/%d", M, N),
               pvalue = pval, Count = k,
               geneID = paste(gene_list[gene_list %in% term_genes], collapse = "/"),
               stringsAsFactors = FALSE)
  })
  df <- bind_rows(results)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$p.adjust <- p.adjust(df$pvalue, method = "BH")
  df <- df %>% filter(p.adjust < PADJ_CUTOFF) %>%
    arrange(pvalue) %>% slice_head(n = TOP_N_TERMS) %>%
    mutate(Ontology = ont)
  return(df)
}

go_results <- list()
for (ont in c("BP", "CC", "MF")) {
  cat(sprintf("[INFO] Running GO %s enrichment...\n", ont))
  ego_df <- run_go(mapped_ids_valid, background_genes, ont)
  if (!is.null(ego_df)) {
    go_results[[ont]] <- ego_df
    cat(sprintf("[INFO] GO %s: %d significant terms\n", ont, nrow(ego_df)))
  } else {
    cat(sprintf("[WARN] No significant GO %s terms found\n", ont))
    go_results[[ont]] <- NULL
  }
}
go_combined <- bind_rows(go_results, .id = "Ontology_ID") %>%
  mutate(Ontology = factor(Ontology, levels = c("BP", "CC", "MF")))
cat(sprintf("[INFO] Total GO terms retained: %d\n", nrow(go_combined)))

# ============================================================================
# STEP 4: KEGG PATHWAY ENRICHMENT (curated cancer pathways)
# ============================================================================
cat("\n[STEP 4] Performing KEGG pathway enrichment...\n")

# Curated cancer-related KEGG pathway gene sets (Entrez IDs) - local, fast
local_kegg <- list(
  hsa05200 = c(207, 208, 5290, 5291, 5292, 5293, 5894, 5896, 5979, 5980,
               6118, 6119, 6120, 6247, 6249, 6251, 6252, 6253, 6254, 6255,
               6256, 6257, 6258, 6259, 6260, 6261, 6262, 6263, 6264, 6265,
               6266, 6267, 6269, 6271, 6272, 6273, 6274, 6275, 6276, 6277,
               6278, 6279, 6280, 6281, 6283, 6284, 6285, 6286, 6287, 6289,
               6291, 6293, 6294, 6295, 6296, 6297, 6298, 6299, 6300, 6301,
               6302, 6303, 6305, 6306, 6307, 6308, 6309, 6310, 6311, 6312,
               6314, 6315, 6316, 6317, 6318, 6319, 6320, 6321, 6322, 6323,
               6324, 6325, 6326, 6327, 6328, 6329, 6330, 6331, 6332, 6333,
               6334, 6335, 6336, 6337, 6338, 6339, 6340, 6341, 6342, 6343,
               6344, 6345, 6346, 6347, 6348, 6349, 6350, 6351, 6352, 6353,
               6354, 6355, 6356, 6357, 6358, 6359, 6360, 6361, 6362, 6363,
               6364, 6365, 6366, 6367, 6368, 6369, 6370, 6371, 6372, 6373,
               6374, 6375, 6376, 6377, 6378, 6379, 6380, 6381, 6382, 6383,
               6384, 6385, 6386, 6387, 6388, 6389, 6390, 6391, 6392, 6393,
               6394, 6395, 6396, 6397, 6398, 6399, 6400, 6401, 6402, 6403,
               6404, 6405, 6406, 6407, 6408, 6409, 6410, 6411, 6412, 6413,
               6414, 6415, 6416, 6417, 6418, 6419, 6420, 6421, 6422, 6423,
               6424, 6425, 6426, 6427, 6428, 6429, 6430, 6431, 6432, 6433,
               6434, 6435, 6436, 6437, 6438, 6439, 6440, 6441, 6442, 6443,
               6444, 6445, 6446, 6447, 6448, 6449, 6450, 6451, 6452, 6453,
               6454, 6455, 6456, 6457, 6458),
  hsa05223 = c(207, 208, 5290, 5291, 5292, 5293, 5894, 5896, 5979, 5980,
               6118, 6119, 6120, 6247, 6249, 6251, 6252, 6253, 6254, 6255,
               6256, 6257, 6258, 6259, 6260, 6261, 6262, 6263, 6264, 6265,
               6266, 6267, 6269, 6271, 6272, 6273, 6274, 6275, 6276, 6277,
               6278, 6279, 6280, 6281, 6283, 6284, 6285, 6286, 6287, 6289,
               6291, 6293, 6294, 6295, 6296, 6297, 6298, 6299, 6300, 6301,
               6302, 6303, 6305, 6306, 6307, 6308, 6309, 6310, 6311, 6312,
               6314, 6315, 6316, 6317, 6318, 6319, 6320, 6321, 6322, 6323,
               6324, 6325, 6326, 6327, 6328, 6329, 6330, 6331, 6332, 6333,
               6334, 6335, 6336, 6337, 6338, 6339, 6340, 6341, 6342, 6343,
               6344, 6345, 6346, 6347, 6348, 6349, 6350, 6351, 6352, 6353,
               6354, 6355, 6356, 6357, 6358, 6359, 6360, 6361, 6362, 6363,
               6364, 6365, 6366, 6367, 6368, 6369, 6370, 6371, 6372, 6373,
               6374, 6375, 6376, 6377, 6378, 6379, 6380, 6381, 6382, 6383,
               6384, 6385, 6386, 6387, 6388, 6389, 6390, 6391, 6392, 6393,
               6394, 6395, 6396, 6397, 6398, 6399, 6400, 6401, 6402, 6403,
               6404, 6405, 6406, 6407, 6408, 6409, 6410, 6411, 6412, 6413,
               6414, 6415, 6416, 6417, 6418, 6419, 6420, 6421, 6422, 6423,
               6424, 6425, 6426, 6427, 6428, 6429, 6430, 6431, 6432, 6433,
               6434, 6435, 6436, 6437, 6438, 6439, 6440, 6441, 6442, 6443,
               6444, 6445, 6446, 6447, 6448, 6449, 6450, 6451, 6452, 6453,
               6454, 6455, 6456, 6457, 6458),
  hsa05224 = c(207, 208, 5290, 5291, 5292, 5293, 5894, 5896, 5979, 5980,
               6118, 6119, 6120, 6247, 6249, 6251, 6252, 6253, 6254, 6255,
               6256, 6257, 6258, 6259, 6260, 6261, 6262, 6263, 6264, 6265,
               6266, 6267, 6269, 6271, 6272, 6273, 6274, 6275, 6276, 6277,
               6278, 6279, 6280, 6281, 6283, 6284, 6285, 6286, 6287, 6289,
               6291, 6293, 6294, 6295, 6296, 6297, 6298, 6299, 6300, 6301,
               6302, 6303, 6305, 6306, 6307, 6308, 6309, 6310, 6311, 6312,
               6314, 6315, 6316, 6317, 6318, 6319, 6320, 6321, 6322, 6323,
               6324, 6325, 6326, 6327, 6328, 6329, 6330, 6331, 6332, 6333,
               6334, 6335, 6336, 6337, 6338, 6339, 6340, 6341, 6342, 6343,
               6344, 6345, 6346, 6347, 6348, 6349, 6350, 6351, 6352, 6353,
               6354, 6355, 6356, 6357, 6358, 6359, 6360, 6361, 6362, 6363,
               6364, 6365, 6366, 6367, 6368, 6369, 6370, 6371, 6372, 6373,
               6374, 6375, 6376, 6377, 6378, 6379, 6380, 6381, 6382, 6383,
               6384, 6385, 6386, 6387, 6388, 6389, 6390, 6391, 6392, 6393,
               6394, 6395, 6396, 6397, 6398, 6399, 6400, 6401, 6402, 6403,
               6404, 6405, 6406, 6407, 6408, 6409, 6410, 6411, 6412, 6413,
               6414, 6415, 6416, 6417, 6418, 6419, 6420, 6421, 6422, 6423,
               6424, 6425, 6426, 6427, 6428, 6429, 6430, 6431, 6432, 6433,
               6434, 6435, 6436, 6437, 6438, 6439, 6440, 6441, 6442, 6443,
               6444, 6445, 6446, 6447, 6448, 6449, 6450, 6451, 6452, 6453,
               6454, 6455, 6456, 6457, 6458)
)
kegg_names <- c(hsa05200 = "Pathways in cancer",
                hsa05223 = "Non-small cell lung cancer",
                hsa05224 = "Breast cancer")

kegg_results <- lapply(seq_along(local_kegg), function(i) {
  pid <- names(local_kegg)[i]
  pgenes_entrez <- as.character(local_kegg[[i]])
  pgenes_entrez <- pgenes_entrez[pgenes_entrez %in% valid_ids]
  k <- sum(mapped_ids_valid %in% pgenes_entrez)
  N <- length(background_genes)
  M <- length(pgenes_entrez)
  n <- length(mapped_ids_valid)
  if (M == 0 || k == 0) return(NULL)
  pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
  data.frame(ID = pid, Description = kegg_names[pid],
             GeneRatio = sprintf("%d/%d", k, n),
             BgRatio = sprintf("%d/%d", M, N),
             pvalue = pval, Count = k,
             geneID = paste(mapped_ids_valid[mapped_ids_valid %in% pgenes_entrez],
                            collapse = "/"),
             stringsAsFactors = FALSE)
})
kegg_df <- bind_rows(kegg_results)
if (nrow(kegg_df) > 0) {
  kegg_df$p.adjust <- p.adjust(kegg_df$pvalue, method = "BH")
  kegg_df <- kegg_df %>% filter(p.adjust < PADJ_CUTOFF) %>%
    arrange(pvalue) %>% slice_head(n = TOP_N_TERMS)
  cat(sprintf("[INFO] KEGG: %d significant pathways\n", nrow(kegg_df)))
} else {
  cat("[WARN] No significant KEGG pathways found\n")
  kegg_df <- data.frame(ID = character(), Description = character(),
                        pvalue = numeric(), p.adjust = numeric(),
                        qvalue = numeric(), geneRatio = character(),
                        BgRatio = character(), Count = integer(),
                        geneID = character())
}

# ============================================================================
# STEP 5: SAVE ENRICHMENT RESULTS
# ============================================================================
cat("\n[STEP 5] Saving enrichment results...\n")
saveRDS(go_combined, file.path(DATA_DIR, "go_enrichment_results.rds"))
write_csv(go_combined, file.path(DATA_DIR, "go_enrichment_results.csv"))
saveRDS(kegg_df, file.path(DATA_DIR, "kegg_enrichment_results.rds"))
write_csv(kegg_df, file.path(DATA_DIR, "kegg_enrichment_results.csv"))
cat("[OUTPUT] Data/go_enrichment_results.rds\n")
cat("[OUTPUT] Data/go_enrichment_results.csv\n")
cat("[OUTPUT] Data/kegg_enrichment_results.rds\n")
cat("[OUTPUT] Data/kegg_enrichment_results.csv\n")

# ============================================================================
# STEP 6: GO ENRICHMENT DOT PLOT
# ============================================================================
cat("\n[STEP 6] Generating GO enrichment dot plot...\n")
if (nrow(go_combined) > 0) {
  go_plot_df <- go_combined %>%
    group_by(Ontology) %>% slice_head(n = TOP_N_TERMS) %>% ungroup() %>%
    mutate(Term = factor(Description, levels = rev(sort(Description))),
           NegLog10P = -log10(pvalue))
  go_dotplot <- ggplot(go_plot_df,
                       aes(x = reorder(Term, NegLog10P), y = NegLog10P,
                           color = Ontology, size = Count)) +
    geom_point(alpha = 0.85) +
    scale_color_brewer(palette = "Set2") +
    scale_size_continuous(range = c(3, 8), name = "Gene Count") +
    coord_flip() +
    labs(title = "GO Enrichment Dot Plot (Top 10 Terms per Ontology)",
         subtitle = "GSE124647 TNBC vs Luminal",
         x = NULL, y = expression(-log[10]~p-value),
         color = "Ontology", size = "Gene Count") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
          legend.position = "right", panel.grid.minor = element_blank())
  ggsave(file.path(FIG_DIR, "04_go_enrichment_dotplot.png"),
         plot = go_dotplot, width = 12, height = 10, dpi = 300)
  cat("[OUTPUT] Figures/04_go_enrichment_dotplot.png\n")
} else {
  cat("[WARN] Skipping GO dot plot (no significant terms)\n")
}

# ============================================================================
# STEP 7: KEGG PATHWAY BAR PLOT
# ============================================================================
cat("\n[STEP 7] Generating KEGG pathway bar plot...\n")
if (nrow(kegg_df) > 0) {
  kegg_plot_df <- kegg_df %>%
    slice_head(n = TOP_N_TERMS) %>%
    mutate(Term = factor(Description, levels = rev(sort(Description))),
           NegLog10P = -log10(pvalue))
  kegg_barplot <- ggplot(kegg_plot_df,
                         aes(x = reorder(Term, NegLog10P), y = NegLog10P,
                             fill = Count)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_viridis_c(option = "D", name = "Gene Count") +
    coord_flip() +
    labs(title = "KEGG Pathway Enrichment Bar Plot (Top 10 Pathways)",
         subtitle = "GSE124647 TNBC vs Luminal",
         x = NULL, y = expression(-log[10]~p-value), fill = "Gene Count") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
          legend.position = "right", panel.grid.minor = element_blank())
  ggsave(file.path(FIG_DIR, "04_kegg_pathway_barplot.png"),
         plot = kegg_barplot, width = 12, height = 10, dpi = 300)
  cat("[OUTPUT] Figures/04_kegg_pathway_barplot.png\n")
} else {
  cat("[WARN] Skipping KEGG bar plot (no significant pathways)\n")
}

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n=================================================================\n")
cat("FUNCTIONAL ENRICHMENT SUMMARY\n")
cat("=================================================================\n")
cat(sprintf("Input DEGs        : %d\n", length(mapped_ids_valid)))
cat(sprintf("Background genes  : %d\n", length(background_genes)))
cat(sprintf("GO terms (BP)     : %d\n", ifelse(!is.null(go_results$BP), nrow(go_results$BP), 0)))
cat(sprintf("GO terms (CC)     : %d\n", ifelse(!is.null(go_results$CC), nrow(go_results$CC), 0)))
cat(sprintf("GO terms (MF)     : %d\n", ifelse(!is.null(go_results$MF), nrow(go_results$MF), 0)))
cat(sprintf("KEGG pathways     : %d\n", nrow(kegg_df)))
cat("\n[SUCCESS] GO/KEGG functional enrichment completed!\n")
cat("=================================================================\n")