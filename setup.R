if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

cran_packages <- c("tidyverse", "ggplot2", "pheatmap", "ggrepel", "survival", "survminer", "igraph")
bioc_packages <- c("GEOquery", "DESeq2", "edgeR", "limma", "clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap", "GSVA")

new_cran <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran)) install.packages(new_cran)

new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(new_bioc)) BiocManager::install(new_bioc)

cat("\n[Success] All core Bioinformatic dependencies installed successfully.\n")