# RNA-Seq Data Analysis Pipeline - Setup Script
# This script installs all required dependencies for the pipeline
# Run this script before executing the analysis pipeline

cat("Installing required packages for Bulk RNA-Seq Analysis Pipeline...\n")

# Core dependency management
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Define required packages
required_packages <- c(
  # Core RNA-Seq Analysis
  "DESeq2",
  "edgeR",
  "limma",
  
  # Data manipulation and visualization
  "ggplot2",
  "pheatmap",
  "ggrepel",
  "RColorBrewer",
  "gridExtra",
  "reshape2",
  "dplyr",
  "tidyr",
  
  # Functional enrichment
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot",
  
  # PCA and QC
  "FactoMineR",
  "factoextra",
  
  # Utilities
  "stringr",
  "readr",
  "data.table"
)

# Install CRAN packages
cran_packages <- required_packages[!grepl("^org\\.", required_packages)]
missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]

if (length(missing_cran) > 0) {
  cat("Installing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# Install Bioconductor packages
bioc_packages <- required_packages[grepl("^(DESeq2|edgeR|limma|clusterProfiler|org\\.Hs\\.eg\\.db|AnnotationDbi|enrichplot|FactoMineR|factoextra)$", required_packages)]
missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]

if (length(missing_bioc) > 0) {
  cat("Installing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, update = FALSE)
}

# Install annotation databases
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  cat("Installing org.Hs.eg.db annotation package...\n")
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

cat("\n✓ All dependencies installed successfully!\n")
cat("✓ You can now run the analysis scripts in the Scripts/ directory.\n")

# Verify installations
cat("\nVerifying installations...\n")
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  ✓ %s\n", pkg))
  } else {
    cat(sprintf("  ✗ %s (installation failed)\n", pkg))
  }
}