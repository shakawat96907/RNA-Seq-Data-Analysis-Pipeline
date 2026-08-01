# 🧬 Breast Cancer Bulk RNA-Seq Transcriptomic Analysis & Prognostic Pipeline (GSE124647)

An automated, end-to-end computational pipeline implemented in R (v4.6.0) for publication-ready Differential Expression Modeling, Subtype Profiling (TNBC vs Non-TNBC), Functional Enrichment, Protein-Protein Interaction (PPI) Networks, Survival Biomarker Prognosis, and Immune Microenvironment Deconvolution using high-throughput Bulk RNA-Seq data.

---

## 📌 Analytical Workflow & Core Modules

This computational framework processes raw count matrices from Breast Cancer cohorts (GSE124647) divided into seven analytical pipelines:

### 1. Count Matrix Preprocessing (`01_count_matrix_preprocessing.R`)
* Programmatic data ingestion of GSE124647 (TNBC vs Luminal/HER2 cohorts).
* Low-count expression filtering ($CPM > 1$), library normalization, and Variance Stabilizing Transformations (VST).
* Sample Quality Control (QC) and Principal Component Analysis (PCA).

### 2. DESeq2 Differential Expression Analysis (`02_deseq2_dge_analysis.R`)
* Differential gene expression modeling using `DESeq2`.
* Strict statistical cutoffs: Adjusted p-value ($p_{adj} < 0.05$) and fold change ($\vert{}\log_2 \text{FC}\vert{} \ge 1.5$).

### 3. Volcano & Heatmap Visualization (`03_volcano_and_heatmap_visualization.R`)
* Volcano plots highlighting driver genes with custom statistical thresholds.
* Hierarchical heatmaps (`pheatmap`) with clinical metadata annotations (Subtypes, Receptor Status).

### 4. GO & KEGG Pathway Enrichment (`04_go_kegg_pathway_enrichment.R`)
* Gene Ontology (GO: BP, CC, MF) and KEGG pathway mapping using `clusterProfiler`.
* Gene Set Enrichment Analysis (GSEA) on Hallmark pathways.

### 5. PPI Network & Hub Gene Identification (`05_protein_network_ppi_hub_genes.R`)
* STRING database network integration.
* Topological centrality assessment (Degree, Betweenness) via `igraph` to discover core BRCA hub genes.

### 6. Clinical Prognosis & Survival Analysis (`06_prognosis_survival_analysis.R`)
* Kaplan-Meier survival curves and Cox proportional hazards regression for biomarker validation.

### 7. Tumor Immune Microenvironment Profiling (`07_immune_infiltration_ssgsea.R`)
* Deconvolution of tumor-infiltrating lymphocytes (TILs) and immune cell types across BRCA subtypes via `ssGSEA`/`GSVA`.

---

## 📁 Repository Structure

```text
RNA-Seq-Data-Analysis-Pipeline/
├── Scripts/
│   ├── 01_count_matrix_preprocessing.R
│   ├── 02_deseq2_dge_analysis.R
│   ├── 03_volcano_and_heatmap_visualization.R
│   ├── 04_go_kegg_pathway_enrichment.R
│   ├── 05_protein_network_ppi_hub_genes.R
│   ├── 06_prognosis_survival_analysis.R
│   └── 07_immune_infiltration_ssgsea.R
├── Figures/
│   ├── .gitkeep
│   └── 01_pca_plot.png
├── Data/
│   └── .gitkeep
├── .gitignore
├── setup.R
└── README.md
```

---

## 📜 Declarations & Contact

* **Author:** Md. Shakawat Hossain
* **Affiliation:** Department of Biochemistry and Molecular Biology, Shahjalal University of Science and Technology (SUST), Sylhet-3114, Bangladesh.
* **Correspondence:** shakawathossain96907@gmail.com
* **Personal Website:** [shakawat-hossain.netlify.app](https://shakawat-hossain.netlify.app/)
* **Professional Networks:** [LinkedIn Profile](https://www.linkedin.com/in/md-shakawat-hossain-372143378/) | [ORCID Profile](https://orcid.org/0009-0008-5050-2275)

---

*Maintained by Md. Shakawat Hossain. For collaboration, pipeline optimization, or transcriptomic data analysis inquiries, reach out via E-mail, LinkedIn, or Website.*
