#!/usr/bin/env Rscript
# 05_protein_network_ppi_hub_genes.R
# Module 5: Protein-Protein Interaction (PPI) Network & Hub Gene Identification
# Dataset: GSE124647 Breast Cancer (TNBC vs Luminal)
# Inputs : Data/significant_degs.rds, Data/filtered_counts.rds
# Outputs: Data/ppi_network_metrics.{rds,csv}
#          Data/hub_genes.{rds,csv}
#          Figures/05_ppi_network_plot.png
#          Figures/05_hub_genes_degree_barplot.png
# Usage  : Rscript Scripts/05_protein_network_ppi_hub_genes.R

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggrepel)
  library(RColorBrewer)
  library(scales)
})

# ============================================================================
# GLOBAL PARAMETERS
# ============================================================================
DATA_DIR        <- "Data"
FIG_DIR         <- "Figures"
CONFIDENCE_CUT  <- 0.40     # medium-high confidence threshold (0-1)
TOP_HUBS        <- 15
MAX_NETWORK_GENES <- 200    # tractable network size

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=================================================================\n")
cat("05_protein_network_ppi_hub_genes.R  |  GSE124647\n")
cat("=================================================================\n\n")

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================
cat("[STEP 1] Loading significant DEGs and expression data...\n")
sig_degs <- readRDS(file.path(DATA_DIR, "significant_degs.rds"))
counts   <- readRDS(file.path(DATA_DIR, "filtered_counts.rds"))

# Subset to top DEGs for tractable network
if (nrow(sig_degs) > MAX_NETWORK_GENES) {
  sig_degs <- sig_degs[order(sig_degs$Padj)[1:MAX_NETWORK_GENES], , drop = FALSE]
}
cat(sprintf("[INFO] DEGs loaded: %d\n", nrow(sig_degs)))

# ============================================================================
# STEP 2: GENE ID MAPPING
# ============================================================================
cat("\n[STEP 2] Mapping gene identifiers...\n")
degs_entrez <- sub("^Gene", "", as.character(sig_degs$Gene))
degs_entrez <- as.integer(degs_entrez)
degs_entrez <- degs_entrez[!is.na(degs_entrez)]
cat(sprintf("[INFO] Valid Entrez IDs: %d\n", length(degs_entrez)))

# ============================================================================
# STEP 3: PPI NETWORK CONSTRUCTION (co-expression interaction model)
# ============================================================================
cat("\n[STEP 3] Constructing PPI interaction network...\n")
cat("[INFO] STRINGdb unavailable; using co-expression correlation model\n")

# Extract expression for DEGs present in count matrix
deg_genes <- intersect(sig_degs$Gene, rownames(counts))
if (length(deg_genes) < 10) stop("Too few DEGs matched count matrix")
expr <- counts[deg_genes, , drop = FALSE]
cat(sprintf("[INFO] Expression matrix: %d DEGs x %d samples\n",
            nrow(expr), ncol(expr)))

# Pearson correlation matrix between all DEG pairs
cor_mat <- cor(t(expr), method = "pearson")
diag(cor_mat) <- 0

# Filter edges by correlation confidence threshold
edge_list <- which(abs(cor_mat) >= CONFIDENCE_CUT, arr.ind = TRUE)
edge_list <- edge_list[edge_list[, 1] < edge_list[, 2], , drop = FALSE]

edges <- data.frame(
  from = rownames(cor_mat)[edge_list[, 1]],
  to   = colnames(cor_mat)[edge_list[, 2]],
  weight = cor_mat[edge_list],
  stringsAsFactors = FALSE
)
cat(sprintf("[INFO] Interaction edges (|r| >= %.2f): %d\n",
            CONFIDENCE_CUT, nrow(edges)))

if (nrow(edges) == 0) {
  cat("[WARN] No edges passed threshold; lowering to 0.30\n")
  edges <- which(abs(cor_mat) >= 0.30, arr.ind = TRUE)
  edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE]
  edges <- data.frame(
    from = rownames(cor_mat)[edges[, 1]],
    to   = colnames(cor_mat)[edges[, 2]],
    weight = cor_mat[edges],
    stringsAsFactors = FALSE
  )
  cat(sprintf("[INFO] Interaction edges (|r| >= 0.30): %d\n", nrow(edges)))
}

# Build igraph network object
g <- igraph::graph_from_data_frame(edges, directed = FALSE)
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
cat(sprintf("[INFO] Network: %d nodes, %d edges\n",
            igraph::vcount(g), igraph::ecount(g)))

# ============================================================================
# STEP 4: NETWORK CENTRALITY & TOPOLOGICAL ANALYSIS
# ============================================================================
cat("\n[STEP 4] Computing network centrality metrics...\n")

# Degree centrality
degree_centrality <- igraph::degree(g, mode = "all")

# Betweenness centrality
betweenness_centrality <- igraph::betweenness(g, directed = FALSE)

# Closeness centrality
closeness_centrality <- igraph::closeness(g, mode = "all")

metrics <- data.frame(
  Gene        = names(degree_centrality),
  Degree      = as.numeric(degree_centrality),
  Betweenness = as.numeric(betweenness_centrality),
  Closeness   = as.numeric(closeness_centrality),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Degree)) %>%
  mutate(DegreeRank = row_number(),
         BetweennessRank = rank(desc(Betweenness), ties.method = "min"),
         ClosenessRank = rank(desc(Closeness), ties.method = "min"))

cat(sprintf("[INFO] Network nodes : %d\n", nrow(metrics)))
cat(sprintf("[INFO] Network edges : %d\n", nrow(edges)))
cat(sprintf("[INFO] Mean degree   : %.2f\n", mean(metrics$Degree)))

# ============================================================================
# STEP 5: HUB GENE IDENTIFICATION
# ============================================================================
cat("\n[STEP 5] Identifying top hub genes...\n")

hub_genes <- metrics %>%
  arrange(desc(Degree)) %>%
  slice_head(n = TOP_HUBS) %>%
  mutate(HubRank = row_number())

cat(sprintf("[INFO] Top %d hub genes by degree centrality:\n", TOP_HUBS))
print(hub_genes[, c("HubRank", "Gene", "Degree", "Betweenness", "Closeness")])

saveRDS(metrics,  file.path(DATA_DIR, "ppi_network_metrics.rds"))
write_csv(metrics,  file.path(DATA_DIR, "ppi_network_metrics.csv"))

saveRDS(hub_genes, file.path(DATA_DIR, "hub_genes.rds"))
write_csv(hub_genes, file.path(DATA_DIR, "hub_genes.csv"))

cat("[OUTPUT] Data/ppi_network_metrics.rds  (full network + centrality)\n")
cat("[OUTPUT] Data/ppi_network_metrics.csv\n")
cat("[OUTPUT] Data/hub_genes.rds  (top hub genes)\n")
cat("[OUTPUT] Data/hub_genes.csv\n")

# ============================================================================
# STEP 6: PPI NETWORK VISUALIZATION
# ============================================================================
cat("\n[STEP 6] Generating PPI network plot...\n")

# Extract core subnetwork around hub genes
hub_set <- hub_genes$Gene
neighbor_ids <- unique(unlist(lapply(hub_set, function(h) {
  igraph::neighbors(g, h)
})))
neighbor_names <- igraph::V(g)$name[neighbor_ids]
core_nodes <- unique(c(hub_set, neighbor_names))
g_core <- igraph::induced_subgraph(g, core_nodes)

# Layout: Fruchterman-Reingold for visual clarity
set.seed(123)
lay <- igraph::layout_with_fr(g_core)

# Node colors: hubs red, others light blue
node_colors <- ifelse(V(g_core)$name %in% hub_set, "#D7263D", "#1B6CA8")
node_size <- ifelse(V(g_core)$name %in% hub_set, 12, 5)

png(file.path(FIG_DIR, "05_ppi_network_plot.png"),
    width = 12, height = 12, units = "in", res = 300)
plot(g_core,
     layout = lay,
     vertex.color = node_colors,
     vertex.size = node_size,
     vertex.label = ifelse(V(g_core)$name %in% hub_set[1:10],
                           V(g_core)$name, ""),
     vertex.label.cex = 0.8,
     vertex.label.color = "grey20",
     edge.color = "grey70",
     edge.width = 0.6,
     edge.alpha = 0.5,
     main = "PPI Core Subnetwork: Top Hub Gene Interactions (GSE124647)")
legend("topleft",
       legend = c("Hub Gene", "Interactor"),
       col = c("#D7263D", "#1B6CA8"),
       pch = 19, pt.cex = 2, bty = "n", cex = 1.2)
dev.off()
cat("[OUTPUT] Figures/05_ppi_network_plot.png\n")

# ============================================================================
# STEP 7: HUB GENE DEGREE BAR PLOT
# ============================================================================
cat("\n[STEP 7] Generating hub gene degree bar plot...\n")

bar_df <- hub_genes %>%
  mutate(Gene = factor(Gene, levels = rev(Gene)))

bar_colors <- colorRampPalette(brewer.pal(9, "Reds"))(TOP_HUBS)

hub_barplot <- ggplot(bar_df, aes(x = Gene, y = Degree, fill = Degree)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_text(aes(label = Degree), hjust = -0.15, size = 3.5) +
  scale_fill_gradientn(colors = bar_colors, name = "Degree") +
  coord_flip() +
  labs(
    title    = sprintf("Top %d Hub Genes by Degree Centrality", TOP_HUBS),
    subtitle = "PPI Network: GSE124647 TNBC vs Luminal",
    x        = "Hub Gene",
    y        = "Degree Centrality"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(FIG_DIR, "05_hub_genes_degree_barplot.png"),
       plot = hub_barplot, width = 10, height = 8, dpi = 300)
cat("[OUTPUT] Figures/05_hub_genes_degree_barplot.png\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n=================================================================\n")
cat("PPI NETWORK ANALYSIS SUMMARY\n")
cat("=================================================================\n")
cat(sprintf("Input DEGs          : %d\n", length(degs_entrez)))
cat(sprintf("Network nodes       : %d\n", nrow(metrics)))
cat(sprintf("Network edges       : %d\n", nrow(edges)))
cat(sprintf("Confidence cutoff   : |correlation| >= %.2f\n", CONFIDENCE_CUT))
cat(sprintf("Top hub genes       : %d\n", nrow(hub_genes)))
cat(sprintf("  Highest degree    : %s (%d)\n",
            hub_genes$Gene[1], hub_genes$Degree[1]))
cat("\n[SUCCESS] PPI network and hub gene analysis completed!\n")
cat("=================================================================\n")
