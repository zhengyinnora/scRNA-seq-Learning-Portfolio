# 🧬 Project: Dissecting the Krt8+ Alveolar Differentiation Intermediate (ADI) in Lung Fibrosis

## 1. Project Background
In lung fibrosis, the regeneration of alveolar epithelium is blocked. Healthy stem cells (**AT2 cells**) fail to differentiate into functional AT1 cells and instead become stalled in a pathological intermediate state known as **Krt8+ ADI (Alveolar Differentiation Intermediate)**. This persistence leads to fibrosis rather than repair.

## 2. Analysis Goal
Using single-cell RNA sequencing (scRNA-seq) data from **GSE141259**, this project aims to:
1.  **Identify** the specific subpopulation of Krt8+ ADI cells.
2.  **Characterize** the molecular features of this "stalled" state compared to healthy AT2 cells.
3.  **Validate** key markers of dedifferentiation and cellular senescence.

## 3. Key Findings
Through differential expression analysis (Seurat workflow), we confirmed that Krt8+ ADI cells exhibit:
* **Dedifferentiation:** Loss of AT2 lineage markers (e.g., *Sftpc*).
* **New Identity:** Gain of injury-associated keratin markers (e.g., *Krt8*).
* **Pathological State:** High expression of senescence (*Cdkn1a/p21*) and stress response markers (*S100a6*, *Clu*), indicating a blocked regenerative process.

---
*Analysis performed by Nora using R/Seurat.*


# 🫁 IPF-Translational-Pipeline: From Single-Cell to Spatial Diagnostics

An end-to-end bioinformatics pipeline for Idiopathic Pulmonary Fibrosis (IPF), bridging cross-species single-cell RNA-seq discovery, drug repurposing, machine learning clinical validation, and spatial transcriptomics mapping.

## 🌟 Project Highlights
* **Cross-Species Integration:** Harmonized human and mouse scRNA-seq data to identify conserved pathogenic cell states (e.g., ADI cells).
* **Drug Repurposing Network:** Extracted target interactions bypassing deprecated packages to construct a customized bipartite network (identifying inhibitors like AMY-101).
* **AI Clinical Diagnostics:** Deployed LASSO and Random Forest to compress 200+ targets into a high-precision diagnostic panel achieving an **AUC of 0.871** on an independent clinical cohort (GSE32537).
* **Spatial Niche Mapping:** Engineered a computational simulation framework to map core biomarkers (e.g., ACTA2) and cell types onto 2D tissue coordinates, overcoming Seurat v5 layer fragmentation.

## 🛠️ Pipeline Structure (20-Step Workflow)
The analysis is divided into sequential R scripts, ensuring reproducibility from raw matrices to final publication-ready figures:

* `01` to `15`: Quality Control, Normalization, Clustering, and Trajectory Inference.
* `16_Cross_Species_Integration.R`: Human-Mouse multi-omics anchoring.
* `17_Drug_Repurposing.R`: Targeted therapeutic screening.
* `18_GEO_Data_Prep.R`: Automated clinical metadata recovery and matrix transformation.
* `19_Machine_Learning_Signatures.R`: Diagnostic model training (LASSO + RF).
* `20_Spatial_Mapping_Simulation.R`: Spatial transcriptomics deconvolution and visualization.

## 📊 Key Results
*(Note: You can insert your generated plots here later, like the ROC curve or Spatial Mapping plot)*
1. Core Diagnostic Panel: `ACTA2`, `GPX3`, `CXCL10`.
2. Model Performance: AUC 0.871 (Validated in GSE32537 cohort).
3. Drug Candidates: Successfully mapped inhibitors to fibrotic hubs.

## 💻 Tech Stack
* **Language:** R (Version 4.5.1)
* **Core Libraries:** `Seurat v5`, `glmnet` (LASSO), `randomForest`, `pROC`, `igraph`.
