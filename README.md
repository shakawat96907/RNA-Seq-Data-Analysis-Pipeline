# 🧬 Bulk RNA-Seq Data Analysis & Functional Enrichment Pipeline

An automated, end-to-end computational pipeline implemented in R for robust Differential Gene Expression (DGE) analysis, publication-ready biostatistical visualizations, and comprehensive functional pathway enrichment profiling using high-throughput RNA-Seq transcriptomic data.

---

## 📌 Workflow Overview & Analytical Modules

This computational framework processes filtered RNA-Seq count matrices and performs end-to-end transcriptomic profiling divided into three core analytical pipelines:

### 1. Differential Gene Expression (DGE) Analysis
* **Engine:** Standardized variance estimation and differential gene expression modeling using `DESeq2` and `edgeR`.
* **Preprocessing:** Low-count gene filtering, sample library size normalization, and dispersion estimation.
* **Statistical Thresholds:** Differentially Expressed Genes (DEGs) are categorized using adjusted p-values ($p_{adj} < 0.05$) and fold-change criteria ($\vert{}\log_2 \text{FC}\vert{} \ge 1.0$).

### 2. High-Dimensional Data Visualization
* **Volcano Plots:** Built using `ggplot2` to highlight significantly up- and down-regulated genes with custom statistical cutoff thresholds.
* **Hierarchical Heatmaps:** Rendered via `pheatmap` using Z-score normalized Variance Stabilizing Transformations (VST) to display sample-wise expression clustering.
* **Principal Component Analysis (PCA):** Sample quality control and batch-effect detection across experimental conditions.

### 3. Functional Pathway & Network Enrichment
* **Gene Ontology (GO):** Over-Representation Analysis (ORA) covering Biological Processes (BP), Cellular Components (CC), and Molecular Functions (MF) via `clusterProfiler`.
* **KEGG & Reactome Mapping:** Annotation of perturbed biological pathways and metabolic cascades using `org.Hs.eg.db`.

---

## 📁 Repository Structure

```text
RNA-Seq-Data-Analysis-Pipeline/
├── Scripts/
│   ├── 01_count_matrix_preprocessing.R
│   ├── 02_deseq2_dge_analysis.R
│   ├── 03_volcano_and_heatmap_visualization.R
│   └── 04_go_kegg_pathway_enrichment.R
├── Figures/
│   └── .gitkeep
├── Data/
│   └── .gitkeep
├── .gitignore
├── setup.R
└── README.md
```

---

## 🚀 Getting Started

### Prerequisites
- R (≥ 4.0.0)
- RStudio (recommended)

### Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/RNA-Seq-Data-Analysis-Pipeline.git
cd RNA-Seq-Data-Analysis-Pipeline
```

2. Run the setup script to install dependencies:
```r
source("setup.R")
```

3. Prepare your data:
   - Place your raw count matrix in `Data/` directory
   - Update sample metadata in the appropriate script

---

## 📊 Pipeline Usage

### Step 1: Count Matrix Preprocessing
```r
source("Scripts/01_count_matrix_preprocessing.R")
```
- Filters low-count genes
- Normalizes library sizes
- Generates quality control metrics

### Step 2: DESeq2 Differential Gene Expression Analysis
```r
source("Scripts/02_deseq2_dge_analysis.R")
```
- Performs DGE analysis using DESeq2/edgeR
- Identifies significantly differentially expressed genes
- Outputs results with adjusted p-values and fold changes

### Step 3: Visualization
```r
source("Scripts/03_volcano_and_heatmap_visualization.R")
```
- Generates volcano plots
- Creates hierarchical heatmaps
- Performs PCA analysis

### Step 4: Functional Enrichment Analysis
```r
source("Scripts/04_go_kegg_pathway_enrichment.R")
```
- GO enrichment (BP, CC, MF)
- KEGG pathway analysis
- Generates enrichment plots

---

## 📦 Dependencies

Key R packages used in this pipeline:
- **DESeq2** - Differential expression analysis
- **edgeR** - Empirical analysis of digital gene expression data
- **ggplot2** - Data visualization
- **pheatmap** - Heatmap visualization
- **clusterProfiler** - Functional enrichment analysis
- **org.Hs.eg.db** - Human gene annotation

See `setup.R` for complete dependency list.

---

## 📝 License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

---

## 📧 Contact

For questions or collaborations, please open an issue on GitHub.