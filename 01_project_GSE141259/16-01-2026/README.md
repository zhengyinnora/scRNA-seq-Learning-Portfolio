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
