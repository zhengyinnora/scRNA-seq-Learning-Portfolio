# 🫁 Lung Fibrosis scRNA-seq Analysis (GSE141259) - Day 1

**Date:** 16-1-2026
**Project:** Lung Regeneration & Fibrosis (Bleomycin Model)  
**Author:** Nora

---

## 📝 Summary (摘要)
This document records the progress of Day 1 analysis using the `GSE141259` dataset. The workflow covers raw data preprocessing, quality control, integration of author-provided annotations, and preliminary investigation of cell type dynamics during disease progression (PBS vs. Bleomycin time course).
本文档记录了基于 `GSE141259` 数据集第一天的分析进度。工作流程涵盖了原始数据预处理、质量控制、作者注释信息的整合，以及对疾病进程（PBS 对照组与博来霉素损伤组的时间序列）中细胞类型动态的初步探究。

---

## 📂 Repository Structure (文件结构)
* **`01_preprocessing_standard.R`**: Script for loading raw counts, QC, normalization, and metadata integration. (数据加载、质控、标准化及元数据整合脚本)
* **`02_disease_comparison.R`**: Script for visualizing disease progression and validating marker genes. (疾病进程可视化及标志物验证脚本)

---

## 🛠️ Workflow & Achievements (工作流程与成果)

### Part 1: Data Preprocessing & Authorization (数据预处理与身份确立)
* **Setup Seurat Object**: Successfully loaded the raw count matrix and initialized the Seurat object.
    * *成功加载原始计数矩阵并初始化 Seurat 对象。*
* **Metadata Correction**: Corrected `orig.ident` from the auto-parsed "muc4169" to the project label "lianxixi".
    * *修正元数据：将自动解析的 "muc4169" 标签更改为项目标签 "lianxixi"。*
* **Annotation Integration**: Loaded author-provided metadata (`GSE141259_WholeLung_cellinfo.csv`) and mapped cell types to the object. Handled `NA` values by assigning them as "Unknown".
    * *注释整合：加载作者提供的元数据文件，将细胞类型映射到对象中，并将缺失值标记为 "Unknown"。*
* **Dimensionality Reduction**: Performed standard normalization, scaling, PCA, and UMAP visualization.
    * *降维分析：完成标准化的归一化、缩放、PCA 及 UMAP 可视化。*

### Part 2: Disease Progression Analysis (疾病进程分析)
* **Time-Course Grouping**: Identified 7 experimental groups: `PBS` (Control) and `d3, d7, d10, d14, d21, d28` (Bleomycin injury model).
    * *时间序列分组：识别出 7 个实验组：PBS（对照）及 d3-d28（博来霉素损伤模型）。*
* **Dynamics Visualization**:
    * **Split UMAP**: Visualized cell shifts across time points. (*拆分 UMAP 图：可视化细胞随时间的迁移。*)
    * **Trend Lines**: Plotted proportion changes for key cell types (AM, AT2, Krt8 ADI). (*折线图：绘制关键细胞类型的比例变化趋势。*)
* **Marker Validation (Krt8 Story)**:
    * Verified `Sftpc` (AT2 marker) loss and `Krt8` (injury marker) emergence using **FeaturePlots**.
    * Confirmed `Krt8` expression specifically in the **Krt8 ADI** population (and normal airway cells) but NOT in healthy AT2 cells using **Violin Plots**.
    * *标志物验证：使用特征图验证了 Sftpc 的丢失和 Krt8 的出现；使用小提琴图证实 Krt8 特异性表达于损伤后的 Krt8 ADI 细胞群（及正常气管细胞），而在健康 AT2 细胞中不表达。*

---

## 📊 Key Findings (关键发现)
1.  **Cellular Shift**: A massive infiltration of macrophages and loss of AT2 cells were observed starting from day 3.
    * *细胞演变：从第 3 天开始观察到巨噬细胞的大量浸润和 AT2 细胞的丢失。*
2.  **Emergence of Krt8+ Cells**: A distinct `Krt8 ADI` population appears around d10-d14, coinciding with the peak of tissue injury.
    * *Krt8+ 细胞的出现：在 d10-d14（组织损伤高峰期）观察到明显的 `Krt8 ADI` 细胞群。*
3.  **Molecular Validation**: Confirmed that `Krt8 ADI` cells express high levels of *Krt8*, distinguishing them from healthy alveolar type 2 cells.
    * *分子验证：证实 `Krt8 ADI` 细胞高表达 *Krt8*，将其与健康的肺泡二型细胞区分开来。*

---

## 🔜 Next Steps (下一步计划)
* Perform Differential Expression Analysis (DEG) to identify genes driving the AT2-to-Krt8 transition.
    * *进行差异表达分析 (DEG)，寻找驱动 AT2 向 Krt8 状态转变的关键基因。*

---
# 📂 Work Log: Differential Expression Analysis (Krt8 ADI vs AT2)

**Date:** 2026-01-17  
**Author:** Nora  
**Script:** `03_scripts_runs/03_DEG_analysis.R`  
**Data Source:** GSE141259 (Mouse Lung Fibrosis Model)

---

## 🎯 1. Objective (分析目标)
To identify the molecular mechanisms driving the transition from healthy **AT2 cells** to the pathological **Krt8 ADI (Alveolar Differentiation Intermediate)** state. Specifically, to find Differentially Expressed Genes (DEGs) that characterize the "dedifferentiation" and "stress response" process.

## 🛠️ 2. Workflow (分析流程)
* **Data Loading:** Loaded the processed Seurat object (`lung_obj_final_analysis.rds`).
* **Set Identity:** Switched active identity to `cell.type`.
* **DEG Calculation:** Used `FindMarkers()` to compare `Krt8 ADI` (Group 1) vs `AT2 cells` (Group 2).
* **Visualization:** Generated Volcano Plots and UMAP Feature Plots to validate findings spatially.

## 🧬 3. Key Biological Findings (生物学发现)

Through statistical analysis and spatial visualization, three key dimensions of the Krt8 ADI cell state were identified:

### A. Identity Switch (身份互换)
* **Loss of AT2 Marker:** `Sftpc` (Surfactant Protein C) was significantly **downregulated** (Log2FC ≈ -2.5), indicating a loss of normal alveolar function (Dedifferentiation).
* **Gain of Injury Marker:** `Krt8` (Cytokeratin 8) was significantly **upregulated** (Log2FC ≈ 3.16), serving as the distinct marker for this injury-associated cell state.

### B. Cellular Senescence (细胞衰老/停滞)
* **Cycle Arrest:** `Cdkn1a` (p21) was highly upregulated (Top 10 upregulated genes). This suggests the cells are locked in a senescent state (cell cycle arrest) and unable to complete regeneration.

### C. Stress Response (压力应激)
* **High Stress Level:** `S100a6` and `Clu` (Clusterin) showed the most dramatic upregulation (Log2FC > 5). This reflects an intense survival response to tissue injury.

## 📊 4. Visual Evidence (结果图表)

| Plot Type | Filename | Description |
| :--- | :--- | :--- |
| **Volcano Plot** | `lianxixi_Volcano_Krt8_Final.png` | Highlights the global shift: Sftpc (Left/Down) vs Krt8/S100a6/Cdkn1a (Right/Up). |
| **Feature Plot** | `lianxixi_FeaturePlot_KeyGenes.png` | Spatial confirmation: `Krt8` expression is exclusive to the specific "bridge" population where `Sftpc` expression is lost, co-localizing perfectly with stress markers `S100a6` and `Cdkn1a`. |

## 📝 5. Conclusion (总结)
The analysis confirms that **Krt8 ADI cells** are not merely "different" AT2 cells, but a distinct, pathological cell state characterized by **dedifferentiation (loss of Sftpc)**, **senescence (Cdkn1a)**, and **high stress (S100a6)**. They represent a "stalled" regeneration intermediate in lung fibrosis.

---
*Created with R Seurat v4/v5*

---

# 🧬 Work Log: Functional Enrichment Analysis (GO:BP)

**Date:** 2026-01-19  
**Author:** Nora  
**Script:** `03_scripts_runs/04_Enrichment_analysis.R`  
**Input Data:** Up-regulated DEGs from Krt8 ADI cells (vs AT2)

## 🎯 1. Objective (分析目标)
To decode the functional state of **Krt8 ADI cells** by mapping the previously identified up-regulated genes (e.g., *Krt8*, *S100a6*, *Clu*) to biological pathways using Gene Ontology (GO) enrichment analysis. We aim to answer: "What are these damaged cells actively *doing*?"

## 🛠️ 2. Methodology (方法)
* **Filtering:** Selected significant up-regulated genes (`adj.P.Val < 0.05` & `log2FC > 0.5`) from the DEG table.
* **Annotation:** Converted Gene Symbols to Entrez IDs using `org.Mm.eg.db`.
* **Enrichment:** Performed GO Biological Process (BP) enrichment using `clusterProfiler::enrichGO`.
* **Visualization:** Generated a Dotplot to visualize the top 15 most significant pathways.

## 🧬 3. Key Findings (核心发现)
The analysis reveals that Krt8 ADI cells are metabolically hyper-active despite being cell-cycle arrested.

### A. Hyper-Biosynthesis (疯狂的合成代谢)
* **Top Terms:** `cytoplasmic translation`, `ribosome biogenesis`, `translation at synapse`.
* **Interpretation:** The most significant function is **protein synthesis**. This explains the high expression of structural proteins (*Krt8*) and stress-response proteins (*S100a6*, *Clu*) identified in the previous step. The cells are effectively "factories" running at full capacity to produce survival factors.

### B. High Metabolic Demand (高能量消耗)
* **Key Terms:** `aerobic respiration`, `oxidative phosphorylation`, `ATP synthesis`.
* **Interpretation:** The protein synthesis machinery requires immense energy. The enrichment of respiratory pathways confirms that these cells are burning fuel (ATP) aggressively to maintain their "stalled" but highly active stress state.

## 📊 4. Visual Evidence (结果图表)

| Plot Type | Filename | Description |
| :--- | :--- | :--- |
| **Dotplot** | `GO_Enrichment_Dotplot.png` | Shows the top enriched biological processes. The dominance of "Translation" (Red/Large bubbles) provides strong evidence for the high-biosynthetic state of Krt8 ADI cells. |

## 📝 5. Integrated Conclusion (综合结论)
Combining the DEG results (Jan 17) with today's Enrichment results (Jan 19):
**Krt8 ADI cells** are defined by a paradox:
1.  **Stalled Growth:** They are not dividing (high *Cdkn1a*).
2.  **Hyper-Active Metabolism:** They are actively synthesizing proteins and generating ATP (GO results).
This confirms they are in a **Senescence-Associated Secretory Phenotype (SASP)-like state**, prioritizing survival and stress signaling over regeneration.

---
*Created with clusterProfiler / R*

---
# 🧬 Work Log: Trajectory & Pseudotime Analysis (Monocle3)

**Date:** 2026-01-26  
**Author:** Nora  
**Script:** `05_Trajectory_Analysis.R`  
**Input Data:** Seurat Object (`lung_obj_final_analysis.rds`) - Subset of AT2 & ADI lineages.

---

## 🎯 1. Objective (分析目标)

To reconstruct the continuous developmental trajectory of alveolar regeneration and quantify the transdifferentiation process from healthy **AT2 cells** to the injury-induced **Krt8+ ADI** state.
We aim to answer: *"Is there a continuous path connecting these two states, and what are the driver genes orchestrating this transition?"*

---

## 🛠️ 2. Methodology (方法与修正)

The analysis was performed using **Monocle3** with the following key steps:

* **Data Structure Conversion:** Successfully converted the Seurat object to Monocle3 `cell_data_set` (CDS).
    * *Correction:* Manually transferred `cell_type` metadata to resolve the ID matching error.
* **Trajectory Inference:** Learned the principal graph on the UMAP embedding.
    * *Result:* A continuous "main path" was identified connecting the AT2 cluster to the ADI cluster.
* **Pseudotime Calculation:**
    * Defined the **Root Node** based on the highest density of `AT2 cells`.
    * Calculated pseudotime values for all cells (0 = Start/AT2, >50 = End/ADI).
* **Driver Gene Identification:**
    * Used **Moran’s I test** (`graph_test`) to find genes with spatially coherent expression along the trajectory.
    * Grouped these genes into **Co-expression Modules** to identify the ADI-specific gene signature.

---

## 🧬 3. Key Findings (核心发现)

### A. Confirmation of Lineage Plasticity (细胞可塑性验证)
The trajectory analysis confirms a direct, continuous lineage relationship between AT2 cells and Krt8 ADI cells. The process is not discrete but shows a gradual transition, supporting the hypothesis of **AT2-to-ADI transdifferentiation**.

### B. Pseudotime Gradient (拟时序梯度)
* **Early Phase:** Corresponds to naive AT2 cells (High expression of surfactant genes).
* **Late Phase:** Corresponds to the Krt8 ADI state, characterized by the loss of AT2 identity and acquisition of a stress/injury phenotype.

### C. Identification of the "ADI Module" (锁定关键基因团伙)
We successfully identified a specific gene module that is sharply upregulated at the trajectory terminus (ADI state).
* **Key Markers identified:** `Krt8`, `Krt19`, `Lgals3`.
* **Biological Implication:** This module represents the core transcriptional machinery driving the injury response and preventing re-differentiation (the "stalled" state).

---

## 📊 4. Visual Evidence (结果图表)

| Plot Type | Filename | Description |
| :--- | :--- | :--- |
| **Trajectory** | `Monocle3_Trajectory_Initial.png` | Visualization of the learned black line (principal graph) on UMAP. |
| **Pseudotime** | `Monocle3_Pseudotime_Final.png` | Cells colored by pseudotime (Purple -> Yellow) showing the direction of differentiation. |
| **Heatmap** | `Monocle3_Heatmap_Modules.png` | Heatmap of gene modules, revealing the specific group of genes activated at the end of the trajectory. |
| **Gene List** | `Krt8_ADI_Module_Genes.csv` | **[Important]** The full list of genes co-expressed with Krt8 (The "Suspect List"). |

---

## 📝 Next Steps
* Perform functional enrichment analysis (GO/KEGG) on the identified "Krt8 Module" to understand the biological functions (e.g., wound healing, senescence).

---
# 📅 2026-01-27: Functional Enrichment Analysis (GO Terms)

**Script:** `06_Module_Enrichment.R`
**Input Data:** `Krt8_ADI_Module_Genes.csv` (The target gene list identified in Step 05)
**Status:** ✅ Completed

## 🎯 1. Objective (分析目标)
To decode the biological function of the "Krt8 ADI Module".
We know *which* genes are in this module (e.g., *Krt8*, *Lgals3*), but we need to understand *what* cellular processes they orchestrate.
* **Question:** Are these cells dying? Proliferating? Or differentiating?

## 🛠️ 2. Methodology (方法)
* **Tool:** `clusterProfiler` (R package).
* **Database:** `org.Mm.eg.db` (Mouse Genome).
* **Analysis Type:** Gene Ontology (GO) Enrichment - **Biological Process (BP)**.
* **Process:**
    1.  Converted Gene Symbols to Entrez IDs.
    2.  Performed enrichment test (p-value cutoff < 0.05).
    3.  Visualized top 15 enriched terms using a **Dotplot**.

## 🧬 3. Key Findings (核心发现)
The analysis reveals a distinct "Reprogramming & Repair" signature:

1.  **Squamous Metaplasia (鳞状化生):**
    * **Top Terms:** `epidermis development`, `skin development`.
    * **Interpretation:** The lung AT2 cells are adopting "skin-like" properties (expression of Keratins). This is a classic stress response to thicken the alveolar barrier and protect against injury (forming a "callus").
2.  **Barrier Reconstruction (屏障重建):**
    * **Top Terms:** `tight junction assembly`, `cell-cell junction organization`.
    * **Interpretation:** The cells are actively building connections to seal the injured alveoli and prevent leakage.
3.  **Regenerative Potential (再生潜力):**
    * **Top Terms:** `lung morphogenesis`, `epithelial tube branching`.
    * **Interpretation:** Despite the injury phenotype, the cells retain the developmental memory required to regenerate lung structure.

## 📊 4. Output Files (产出)
* 🖼️ **Visualization:** `Krt8_Module_GO_Enrichment.png` (The Dotplot).
* 📋 **Full Table:** `Krt8_Module_GO_Table.csv` (Detailed statistics for all terms).

> **Conclusion:** The Krt8 ADI state represents a **"Barrier-Strengthening & Repair"** phase, characterized by the transient activation of epidermal programs to survive acute injury.
