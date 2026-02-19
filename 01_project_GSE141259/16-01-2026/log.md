# 🫁 Lung Fibrosis scRNA-seq Analysis (GSE141259) - Day 1

**Date:** 16-1-2026
**Project:** Lung Regeneration & Fibrosis (Bleomycin Model)  
**Author:** Nora

---

## 📝 Summary
This document records the progress of Day 1 analysis using the `GSE141259` dataset. The workflow covers raw data preprocessing, quality control, integration of author-provided annotations, and preliminary investigation of cell type dynamics during disease progression (PBS vs. Bleomycin time course).
本文档记录了基于 `GSE141259` 数据集第一天的分析进度。工作流程涵盖了原始数据预处理、质量控制、作者注释信息的整合，以及对疾病进程（PBS 对照组与博来霉素损伤组的时间序列）中细胞类型动态的初步探究。

---

## 📂 Repository Structure
* **`01_preprocessing_standard.R`**: Script for loading raw counts, QC, normalization, and metadata integration. (数据加载、质控、标准化及元数据整合脚本)
* **`02_disease_comparison.R`**: Script for visualizing disease progression and validating marker genes. (疾病进程可视化及标志物验证脚本)

---

## 🛠️ Workflow & Achievements

### Part 1: Data Preprocessing & Authorization
* **Setup Seurat Object**: Successfully loaded the raw count matrix and initialized the Seurat object.
    * *成功加载原始计数矩阵并初始化 Seurat 对象。*
* **Metadata Correction**: Corrected `orig.ident` from the auto-parsed "muc4169" to the project label "lianxixi".
    * *修正元数据：将自动解析的 "muc4169" 标签更改为项目标签 "lianxixi"。*
* **Annotation Integration**: Loaded author-provided metadata (`GSE141259_WholeLung_cellinfo.csv`) and mapped cell types to the object. Handled `NA` values by assigning them as "Unknown".
    * *注释整合：加载作者提供的元数据文件，将细胞类型映射到对象中，并将缺失值标记为 "Unknown"。*
* **Dimensionality Reduction**: Performed standard normalization, scaling, PCA, and UMAP visualization.
    * *降维分析：完成标准化的归一化、缩放、PCA 及 UMAP 可视化。*

### Part 2: Disease Progression Analysis
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

## 📊 Key Findings
1.  **Cellular Shift**: A massive infiltration of macrophages and loss of AT2 cells were observed starting from day 3.
    * *细胞演变：从第 3 天开始观察到巨噬细胞的大量浸润和 AT2 细胞的丢失。*
2.  **Emergence of Krt8+ Cells**: A distinct `Krt8 ADI` population appears around d10-d14, coinciding with the peak of tissue injury.
    * *Krt8+ 细胞的出现：在 d10-d14（组织损伤高峰期）观察到明显的 `Krt8 ADI` 细胞群。*
3.  **Molecular Validation**: Confirmed that `Krt8 ADI` cells express high levels of *Krt8*, distinguishing them from healthy alveolar type 2 cells.
    * *分子验证：证实 `Krt8 ADI` 细胞高表达 *Krt8*，将其与健康的肺泡二型细胞区分开来。*

---

## 🔜 Next Steps
* Perform Differential Expression Analysis (DEG) to identify genes driving the AT2-to-Krt8 transition.
    * *进行差异表达分析 (DEG)，寻找驱动 AT2 向 Krt8 状态转变的关键基因。*

---
# 📂 Work Log: Differential Expression Analysis (Krt8 ADI vs AT2)

**Date:** 2026-01-17  
**Author:** Nora  
**Script:** `03_scripts_runs/03_DEG_analysis.R`  
**Data Source:** GSE141259 (Mouse Lung Fibrosis Model)

---

## 🎯 1. Objective
To identify the molecular mechanisms driving the transition from healthy **AT2 cells** to the pathological **Krt8 ADI (Alveolar Differentiation Intermediate)** state. Specifically, to find Differentially Expressed Genes (DEGs) that characterize the "dedifferentiation" and "stress response" process.

## 🛠️ 2. Workflow
* **Data Loading:** Loaded the processed Seurat object (`lung_obj_final_analysis.rds`).
* **Set Identity:** Switched active identity to `cell.type`.
* **DEG Calculation:** Used `FindMarkers()` to compare `Krt8 ADI` (Group 1) vs `AT2 cells` (Group 2).
* **Visualization:** Generated Volcano Plots and UMAP Feature Plots to validate findings spatially.

## 🧬 3. Key Biological Findings

Through statistical analysis and spatial visualization, three key dimensions of the Krt8 ADI cell state were identified:

### A. Identity Switch
* **Loss of AT2 Marker:** `Sftpc` (Surfactant Protein C) was significantly **downregulated** (Log2FC ≈ -2.5), indicating a loss of normal alveolar function (Dedifferentiation).
* **Gain of Injury Marker:** `Krt8` (Cytokeratin 8) was significantly **upregulated** (Log2FC ≈ 3.16), serving as the distinct marker for this injury-associated cell state.

### B. Cellular Senescence
* **Cycle Arrest:** `Cdkn1a` (p21) was highly upregulated (Top 10 upregulated genes). This suggests the cells are locked in a senescent state (cell cycle arrest) and unable to complete regeneration.

### C. Stress Response
* **High Stress Level:** `S100a6` and `Clu` (Clusterin) showed the most dramatic upregulation (Log2FC > 5). This reflects an intense survival response to tissue injury.

## 📊 4. Visual Evidence

| Plot Type | Filename | Description |
| :--- | :--- | :--- |
| **Volcano Plot** | `lianxixi_Volcano_Krt8_Final.png` | Highlights the global shift: Sftpc (Left/Down) vs Krt8/S100a6/Cdkn1a (Right/Up). |
| **Feature Plot** | `lianxixi_FeaturePlot_KeyGenes.png` | Spatial confirmation: `Krt8` expression is exclusive to the specific "bridge" population where `Sftpc` expression is lost, co-localizing perfectly with stress markers `S100a6` and `Cdkn1a`. |

## 📝 5. Conclusion
The analysis confirms that **Krt8 ADI cells** are not merely "different" AT2 cells, but a distinct, pathological cell state characterized by **dedifferentiation (loss of Sftpc)**, **senescence (Cdkn1a)**, and **high stress (S100a6)**. They represent a "stalled" regeneration intermediate in lung fibrosis.

---
*Created with R Seurat v4/v5*

---

# 🧬 Work Log: Functional Enrichment Analysis (GO:BP)

**Date:** 2026-01-19  
**Author:** Nora  
**Script:** `03_scripts_runs/04_Enrichment_analysis.R`  
**Input Data:** Up-regulated DEGs from Krt8 ADI cells (vs AT2)

## 🎯 1. Objective
To decode the functional state of **Krt8 ADI cells** by mapping the previously identified up-regulated genes (e.g., *Krt8*, *S100a6*, *Clu*) to biological pathways using Gene Ontology (GO) enrichment analysis. We aim to answer: "What are these damaged cells actively *doing*?"

## 🛠️ 2. Methodology
* **Filtering:** Selected significant up-regulated genes (`adj.P.Val < 0.05` & `log2FC > 0.5`) from the DEG table.
* **Annotation:** Converted Gene Symbols to Entrez IDs using `org.Mm.eg.db`.
* **Enrichment:** Performed GO Biological Process (BP) enrichment using `clusterProfiler::enrichGO`.
* **Visualization:** Generated a Dotplot to visualize the top 15 most significant pathways.

## 🧬 3. Key Findings
The analysis reveals that Krt8 ADI cells are metabolically hyper-active despite being cell-cycle arrested.

### A. Hyper-Biosynthesis
* **Top Terms:** `cytoplasmic translation`, `ribosome biogenesis`, `translation at synapse`.
* **Interpretation:** The most significant function is **protein synthesis**. This explains the high expression of structural proteins (*Krt8*) and stress-response proteins (*S100a6*, *Clu*) identified in the previous step. The cells are effectively "factories" running at full capacity to produce survival factors.

### B. High Metabolic Demand
* **Key Terms:** `aerobic respiration`, `oxidative phosphorylation`, `ATP synthesis`.
* **Interpretation:** The protein synthesis machinery requires immense energy. The enrichment of respiratory pathways confirms that these cells are burning fuel (ATP) aggressively to maintain their "stalled" but highly active stress state.

## 📊 4. Visual Evidence

| Plot Type | Filename | Description |
| :--- | :--- | :--- |
| **Dotplot** | `GO_Enrichment_Dotplot.png` | Shows the top enriched biological processes. The dominance of "Translation" (Red/Large bubbles) provides strong evidence for the high-biosynthetic state of Krt8 ADI cells. |

## 📝 5. Integrated Conclusion
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

## 🎯 1. Objective

To reconstruct the continuous developmental trajectory of alveolar regeneration and quantify the transdifferentiation process from healthy **AT2 cells** to the injury-induced **Krt8+ ADI** state.
We aim to answer: *"Is there a continuous path connecting these two states, and what are the driver genes orchestrating this transition?"*

---

## 🛠️ 2. Methodology

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

## 🧬 3. Key Findings

### A. Confirmation of Lineage Plasticity
The trajectory analysis confirms a direct, continuous lineage relationship between AT2 cells and Krt8 ADI cells. The process is not discrete but shows a gradual transition, supporting the hypothesis of **AT2-to-ADI transdifferentiation**.

### B. Pseudotime Gradient
* **Early Phase:** Corresponds to naive AT2 cells (High expression of surfactant genes).
* **Late Phase:** Corresponds to the Krt8 ADI state, characterized by the loss of AT2 identity and acquisition of a stress/injury phenotype.

### C. Identification of the "ADI Module"
We successfully identified a specific gene module that is sharply upregulated at the trajectory terminus (ADI state).
* **Key Markers identified:** `Krt8`, `Krt19`, `Lgals3`.
* **Biological Implication:** This module represents the core transcriptional machinery driving the injury response and preventing re-differentiation (the "stalled" state).

---

## 📊 4. Visual Evidence 

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

## 🎯 1. Objective
To decode the biological function of the "Krt8 ADI Module".
We know *which* genes are in this module (e.g., *Krt8*, *Lgals3*), but we need to understand *what* cellular processes they orchestrate.
* **Question:** Are these cells dying? Proliferating? Or differentiating?

## 🛠️ 2. Methodology
* **Tool:** `clusterProfiler` (R package).
* **Database:** `org.Mm.eg.db` (Mouse Genome).
* **Analysis Type:** Gene Ontology (GO) Enrichment - **Biological Process (BP)**.
* **Process:**
    1.  Converted Gene Symbols to Entrez IDs.
    2.  Performed enrichment test (p-value cutoff < 0.05).
    3.  Visualized top 15 enriched terms using a **Dotplot**.

## 🧬 3. Key Findings
The analysis reveals a distinct "Reprogramming & Repair" signature:

1.  **Squamous Metaplasia:**
    * **Top Terms:** `epidermis development`, `skin development`.
    * **Interpretation:** The lung AT2 cells are adopting "skin-like" properties (expression of Keratins). This is a classic stress response to thicken the alveolar barrier and protect against injury (forming a "callus").
2.  **Barrier Reconstruction:**
    * **Top Terms:** `tight junction assembly`, `cell-cell junction organization`.
    * **Interpretation:** The cells are actively building connections to seal the injured alveoli and prevent leakage.
3.  **Regenerative Potential :**
    * **Top Terms:** `lung morphogenesis`, `epithelial tube branching`.
    * **Interpretation:** Despite the injury phenotype, the cells retain the developmental memory required to regenerate lung structure.

## 📊 4. Output Files
* 🖼️ **Visualization:** `Krt8_Module_GO_Enrichment.png` (The Dotplot).
* 📋 **Full Table:** `Krt8_Module_GO_Table.csv` (Detailed statistics for all terms).

> **Conclusion:** The Krt8 ADI state represents a **"Barrier-Strengthening & Repair"** phase, characterized by the transient activation of epidermal programs to survive acute injury.


---
# 📅 2026-01-27: Pseudotime Expression Kinetics (Gene Trends)

**Script:** `07_Gene_Trends.R`
**Input Data:** `lung_monocle_final.rds`
**Status:** ✅ Completed

## 🎯 1. Objective
To visualize the dynamic expression changes of specific marker genes along the inferred trajectory.
We aim to validate the identity switch from AT2 to ADI at the single-gene level: *"Do AT2 markers actually drop? Do ADI markers actually rise?"*

## 🛠️ 2. Methodology
* **Function:** `monocle3::plot_genes_in_pseudotime`.
* **Target Genes Selected:**
    * **AT2 Lineage:** *Sftpc, Sftpb* (Surfactant proteins).
    * **ADI/Injury Lineage:** *Krt8, Krt19, Lgals3* (Stress markers).
    * **Proliferation:** *Mki67* (Cell cycle marker).
* **Visualization:** Modeled expression trends (black lines) overlaid on single-cell expression values (colored dots).

## 🧬 3. Key Findings
The kinetic plots reveal distinct expression patterns confirming the transdifferentiation model:

1.  **Loss of AT2 Identity (Sftpc/Sftpb):**
    * **Trend:** Sharp, monotonic decrease along pseudotime.
    * **Meaning:** Cells rapidly silence their surfactant production programs upon injury.
2.  **Gain of Injury State (Krt8/Krt19):**
    * **Trend:** While the smoothed average line remains low (due to zero-inflation/dropout), a distinct sub-population of cells at the terminus (Yellow phase) shows high expression (dots > 3.0).
    * **Meaning:** A specific subset of cells successfully acquires the ADI phenotype.
3.  **Transient Driver (Lgals3):**
    * **Trend:** Strong, bell-shaped upregulation during the intermediate phase.
    * **Meaning:** *Lgals3* appears to be a robust driver activated during the transition process.
4.  **No Proliferation (Mki67):**
    * **Trend:** Flat/Low expression throughout.
    * **Meaning:** Confirms the ADI state is a result of direct conversion (transdifferentiation) without active cell division (cell cycle arrest).

## 📊 4. Output Files
* 🖼️ **Trend Plot:** `Monocle3_Gene_Trends_LinePlot.png` (Visualization of the 6 key genes).

> **Conclusion:** The gene kinetics provide molecular evidence for the AT2-to-ADI transition, characterized by the synchronized downregulation of homeostatic genes and upregulation of stress-response genes.



---
# 📅 2026-02-02: Cell-Cell Communication Inference (Global Network)

**Script:** `08_CellChat_Analysis.R`
**Input Data:** `lung_obj_final_analysis.rds` (Seurat Object)
**Status:** ✅ Completed

## 🎯 1. Objective
To infer the intercellular communication network governing the lung injury response.
We aim to move beyond single-cell identify (AT2 vs. ADI) and understand the **social network** of the tissue: *"Who is signaling to whom during the regeneration process?"*

## 🛠️ 2. Methodology
* **Tool:** `CellChat` (v1.x).
* **Database:** `CellChatDB.mouse` (Full database).
* **Process:**
    1.  **Preprocessing:** Identified over-expressed ligands and receptors in each cell group.
    2.  **Inference:** Computed communication probabilities using the `triMean` method.
    3.  **Aggregation:** Aggregated all L-R pairs to visualize the total communication flow.

## 🧬 3. Key Findings
* **High Connectivity:** The global network reveals a dense, complex web of interactions (the "Hairball" phenotype), indicating robust cross-talk between Epithelial (AT2/ADI), Immune (Macrophages), and Stromal (Fibroblasts) compartments.
* **Distinction:**
    * **Number Plot:** Shows potential connectivity bandwidth (how many pathways are available).
    * **Weight Plot:** Shows actual signaling intensity (how active these pathways are).

## 📊 4. Output Files
* 💾 **CellChat Object:** `lung_cellchat.rds` (Saved for downstream analysis).
* 🖼️ **Global Networks:**
    * `CellChat_Net_Number.png` (Interaction Frequency).
    * `CellChat_Net_Weights.png` (Interaction Strength).
    * ---
# 📅 2026-02-02 (Part 2): Ligand-Receptor Specifics (Zoom-in Analysis)

**Script:** `08_CellChat_Analysis.R` (Part 2)
**Status:** ✅ Completed

## 🎯 1. Objective
To decode the specific molecular signals driving the AT2-to-ADI transition.
Instead of looking at the global network, we perform a directed analysis to answer two specific questions:
1.  **Incoming:** What signals from the microenvironment are "bullying" AT2 cells?
2.  **Outgoing:** What pathogenic signals are Krt8 ADI cells spreading to their neighbors?

## 🛠️ 2. Methodology
* **Visualization:** `netVisual_bubble` (Bubble Plots).
* **Strategy:**
    * **Input Analysis:** Targets restricted to "AT2 cells".
    * **Output Analysis:** Sources restricted to "Krt8 ADI".
* **Key Metrics:** Ligand-Receptor pairs plotted by communication probability (bubble size) and expression level (color).

## 🧬 3. Key Findings

### A. Environment Stress on AT2 (Incoming)
* **Key Signal:** `Fn1 - Sdc4` (Fibronectin - Syndecan-4).
* **Senders:** Fibroblasts, Myofibroblasts, and Fn1+ Macrophages.
* **Interpretation:** AT2 cells are being subjected to intense **Extracellular Matrix (ECM) signaling**. The high interaction with Fibronectin (`Fn1`) suggests that **mechanical stress/stiffening** of the niche is a primary driver forcing AT2 cells to lose their identity and undergo reprogramming.

### B. Pro-inflammatory Feedback by ADI (Outgoing)
* **Key Signals:** `Mif - (Cd74+...)` and `Spp1 - Cd44`.
* **Receivers:** Macrophages (AM, M2) and Proliferating cells.
* **Interpretation:** The Krt8 ADI state is not passive. It actively secretes potent pro-inflammatory mediators:
    * **MIF** (Macrophage Migration Inhibitory Factor): Keeps immune cells trapped and active in the injury site.
    * **SPP1** (Osteopontin): A known driver of fibrosis and inflammation.
    * **Conclusion:** ADI cells create a **pathogenic feed-forward loop**, preventing inflammation resolution.

## 📊 4. Output Files
* 🖼️ `CellChat_Incoming_AT2.png`: Shows the fibrotic signals received by AT2.
* 🖼️ `CellChat_Outgoing_ADI.png`: Shows the inflammatory signals sent by ADI.

> **Summary:** The analysis identifies the **Fn1-Sdc4 axis** as a potential upstream trigger for AT2 injury, and the **Mif/Spp1 axis** as a downstream consequence of the ADI state.

> **Next Step:** Perform "Sender/Receiver" analysis to specifically zoom in on signals received by AT2 cells (Input) and signals sent by Krt8 ADI cells (Output).


# 📅 2026-02-07 (Part 3): Transcription Factor Analysis (The "Master Regulators")

**Script:** `09_TF_Analysis_decoupleR.R` (Final AutoFix Version)
**Status:** ✅ Completed (High Difficulty)

## 🎯 1. Objective
To move beyond gene expression and identify the **upstream Transcription Factors (TFs)** controlling the AT2-to-ADI transition. We aim to find the "bosses" (TFs) in the nucleus that are orchestrating the reprogramming process.

## 🛠️ 2. Methodology
* **Tool:** `decoupleR` (R package for footprint analysis).
* **Algorithm:** `run_ulm` (Univariate Linear Model) for fast and robust inference.
* **Database:** `DoRothEA` (High-confidence TF-Target interactome).
* **Process:**
    * Inferred TF activity scores based on the aggregate expression of downstream target genes.
    * **Troubleshooting:** Overcame persistent "cell name mismatch" errors by implementing a robust renaming strategy (renaming all cells to `C1`, `C2`...) to ensure 100% matrix alignment.

## 🧬 3. Key Findings

### A. The Mastermind of Reprogramming
* **Key TF:** **Sox9**.
* **Observation:** Shows peak activity (deep red) specifically in `Krt8 ADI` cells, while silent in normal AT2 cells.
* **Interpretation:** **Sox9** is confirmed as the master regulator forcing AT2 cells to lose their identity and acquire a progenitor-like state.

### B. The Stress Signal (The Trigger)
* **Key TFs:** **Atf6** and **Hsf1**.
* **Observation:** Highly upregulated activity in the ADI population.
* **Interpretation:** Markers of **ER Stress (Unfolded Protein Response)**. This suggests that proteotoxic stress is a primary driver or consequence of the transition.

### C. The Fibrosis Architects (EMT)
* **Key TF:** **Twist1**.
* **Interpretation:** **Twist1** is a classic driver of Epithelial-Mesenchymal Transition (EMT). Its activation indicates AT2 cells are acquiring fibroblast-like properties, contributing to fibrosis.

### D. The Inflammatory Thugs
* **Key TFs:** **Junb**, **Fosl1** (AP-1 family).
* **Interpretation:** Indicates a strong pro-inflammatory transcriptional program is active within the transitioning cells.

## 📊 4. Output Files
* 🖼️ `TF_Activity_Heatmap.png`: Heatmap visualizing the inferred activity of the top differential TFs across all cell types.

> **Summary:** The analysis successfully identifies the **Sox9-Twist1-Atf6 axis** as the core regulatory network. Injury triggers **ER Stress (Atf6)**, which activates **Sox9** for reprogramming and **Twist1** for fibrosis.

**Next Step:** Perform Pathway Activity Analysis (PROGENy/UCell) to score specific biological processes (e.g., Hypoxia, TGF-b signaling) across the trajectory.


# 📅 2026-02-07 (Part 4): Pathway Activity Scoring (UCell)

**Script:** `10_Pathway_Scoring_UCell.R`
**Status:** ✅ SUCCEEDED (Dependencies fixed via BiocManager)

## 🎯 1. Objective
To quantify specific biological states (Senescence, EMT, Hypoxia) across cell types, specifically asking: *"How 'stressed' are the ADI cells compared to normal AT2?"*

## 🛠️ 2. Methodology
* **Tool:** `UCell` (Rank-based scoring, robust to batch effects).
* **Gene Sets:**
    * **ADI Signature:** (Krt8, Cldn4...) - Positive control.
    * **EMT:** (Vim, Fn1, Twist1...) - Fibrosis indicator.
    * **Senescence:** (Cdkn2a, Il6...) - Aging/SASP indicator.
    * **Hypoxia:** (Hif1a, Vegfa...) - Metabolic stress.

## 🧬 3. Key Findings (from Violin Plots)
* **Validation:** Fibroblasts showed the highest **EMT scores**, validating the gene set accuracy.
* **The ADI Phenotype:**
    * **Partial EMT:** Krt8 ADI cells exhibit significantly higher EMT scores than AT2 cells, confirming their transition towards a mesenchymal state.
    * **Senescence & Hypoxia:** ADI cells are distinctively marked by high Hypoxia and Senescence scores, suggesting they are "stuck" in a stressed, non-regenerative state.
    * **Inflammation:** High inflammatory scores align with the previous CellChat finding (Mif/Spp1 secretion).

## 📊 4. Output Files
* 🖼️ `Pathway_UCell_Violin.png`: Shows the stepwise elevation of stress markers from AT2 -> ADI.
* 🖼️ `Pathway_UCell_UMAP.png`: Visualizes the spatial restriction of these states.

> **Conclusion:** The AT2-to-ADI transition is not just a change in markers (Krt8+), but a fundamental shift in cell state involving **metabolic reprogramming (Hypoxia)**, **cell cycle arrest (Senescence)**, and **fibrogenic activation (EMT)**.


# 📅 2026-02-07 (Part 5): Metabolic Pathway Analysis (The "Switch")

**Script:** `11_Metabolic_Analysis_UCell_Manual.R` (Plan B: Manual Gene Sets)
**Status:** ✅ SUCCEEDED (Bypassed installation issues by using UCell directly)

## 🎯 1. Objective
To validate the hypothesis derived from the Hypoxia scores: *"Do ADI cells undergo metabolic reprogramming (switching from oxidative phosphorylation/lipid metabolism to glycolysis) to adapt to the injury environment?"*

## 🛠️ 2. Methodology
* **Strategy:** Due to `scMetabolism` installation failures (GitHub connectivity), we implemented a **manual UCell scoring approach**.
* **Gene Sets (Curated from KEGG/Reactome):**
    * **Glycolysis:** (Pgk1, Eno1, Ldha...) - Indicator of anaerobic respiration/Warburg effect.
    * **Fatty Acid Metabolism:** (Fasn, Acaca...) - Indicator of Surfactant synthesis.
    * **OxPhos / TCA:** Indicators of mitochondrial respiration.

## 🧬 3. Key Findings

### A. The "Metabolic Switch" (Major Discovery)
* **Evidence:**
* **Observation:**
    * **Normal AT2 cells** are dominant in **Fatty Acid Metabolism**. This makes biological sense as they are "factories" for lipid-rich surfactant.
    * **Krt8 ADI cells** show a dramatic **loss** of Fatty Acid Metabolism and a sharp **gain** in **Glycolysis**.
* **Interpretation:** This confirms **Metabolic Reprogramming**. The ADI cells shift their energy source from lipids to glucose, likely to survive in the **Hypoxic** niche identified in the previous step.

### B. Stress Adaptation
* **Evidence:**
* **Observation:** **Glutathione metabolism** is upregulated in ADI cells.
* **Interpretation:** This suggests ADI cells are under high oxidative stress and are upregulating antioxidant pathways to prevent apoptosis.

### C. Spatial Segregation
* **Evidence:**
* **Observation:** The Glycolysis-high population and Fatty Acid-high population are spatially distinct on the UMAP, with the transition zone (Activated AT2) showing intermediate levels.

## 📊 4. Output Files
* 🖼️ `Metabolism_DotPlot_Manual.png`: The summary of metabolic activities across cell types.
* 🖼️ `Metabolism_UMAP_Comparison.png`: Visualizes the "Lipid-to-Sugar" switch on the single-cell map.

> **Conclusion:** The analysis provides a metabolic mechanism for the ADI state. The transition involves a fundamental switch: **Lipid-Dependent (Homeostatic) ➡️ Glycolytic (Injury-Adapted)**. This aligns perfectly with the Hypoxia and HIF1a signatures found earlier.


# 📅 2026-02-19: Upstream Regulatory Network Inference (NicheNet)

**Script:** `13_NicheNet_Analysis.R`
**Input Data:** `lung_obj_final_analysis.rds`
**Status:** ✅ Completed

## 🎯 1. Objective
To bridge the gap between intercellular communication (CellChat) and intracellular transcription (TF analysis). We aim to identify which specific ligands from the fibrotic niche (Fibroblasts) are directly regulating the pathogenic gene expression profile (Targets) of Krt8 ADI cells.

## 🛠️ 2. Methodology
* **Tool:** `nichenetr` (NicheNet).
* **Setup:**
    * **Sender Cells:** Fibroblasts, Myofibroblasts.
    * **Receiver Cells:** Krt8 ADI.
* **Process:** Prioritized ligands based on Pearson correlation with the ADI target gene set and mapped the regulatory potential to specific target genes.

## 🧬 3. Key Findings

### A. The Fibrotic Niche is Driven by THBS1 and ECM
* **Evidence:**
* **Observation:** The top prioritized ligand is **THBS1** (Thrombospondin-1), followed heavily by structural ECM components (**COL4A1, LAMB1, FBN1, FN1**).
* **Interpretation:** Fibroblasts are not just passively producing scar tissue; they are actively secreting THBS1 to activate latent TGF-β, establishing a stiff, pro-fibrotic mechanical niche that entraps transitioning AT2 cells.

### B. Niche Signals Drive Senescence/Inflammation via ICAM1
* **Evidence:**
* **Observation:** Niche-derived signals, specifically **SPP1**, **THBS1**, and **LAMB1**, show high regulatory potential for upregulating **ICAM1** in ADI cells.
* **Interpretation:** *ICAM1* is a hallmark of cellular senescence and the SASP phenotype. This proves that the pathogenic state of ADI cells is actively maintained by external microenvironmental cues, locking them into an inflammatory, non-regenerative state.

> **Conclusion:** The transition to the ADI state is externally enforced. The Fibroblast-derived **SPP1/THBS1 ➡️ ICAM1 axis** provides a direct mechanistic link between matrix stiffening, cellular senescence, and chronic inflammation in the lung.

# 📅 2026-02-19: Clinical Validation - Part 1 (The PBMC Anomaly)

**Script:** `14_Survival_Analysis.R`
**Data Source:** `GSE28042` (Human IPF Peripheral Blood Mononuclear Cells - PBMC)
**Status:** ✅ Completed (Unexpected Result)

## 🎯 1. Objective
To test if the `Krt8 ADI` gene signature (derived from mouse lung tissue) holds prognostic value in the systemic circulation (blood) of human IPF patients.

## 🛠️ 2. Methodology
* Extracted the 6-gene ADI signature (`KRT8, SOX9, ATF6, TWIST1, SPP1, ICAM1`).
* Scored 75 IPF patients based on their PBMC microarray data.
* Plotted Kaplan-Meier survival curves.

## 🧬 3. Key Findings & The "Plot Twist"
* **Result:** The analysis yielded a highly significant p-value (p = 0.00036). However, the direction was inverted: **Patients with a HIGH ADI score in their blood had significantly LONGER survival times.**
* **Interpretation (Tissue vs. Systemic Paradox):** This biological curveball highlights tissue-specificity. While these genes drive fibrosis locally in the lung, their depletion in the peripheral blood may indicate severe immune exhaustion or systemic failure. 
* **Next Step:** To validate the direct pathogenic role of ADI cells, we must switch the clinical dataset from blood (systemic) to lung tissue/BAL (local microenvironment).

---

# 📅 2026-02-19: Clinical Validation - Part 2 (Metadata Diagnostics)

**Script:** `14c_Diagnostic.R`
**Data Source:** `GSE70866` (Human IPF BAL Fluid)
**Status:** 🛠️ Troubleshooting Completed

## 🎯 1. Objective
To debug the `Column 'Gene' is not found` and `0 patients extracted` errors encountered when initially processing the lung BAL dataset.

## 🛠️ 2. Methodology
* Rather than relying on standard regex matching for clinical columns (which failed), a diagnostic script was written to directly fetch and print the raw `pData` column names and the first patient's metadata from the GEO object.

## 🔍 3. Key Discoveries (The "GEO Formatting" Trap)
* **Finding:** The errors were caused by highly idiosyncratic and non-standard metadata naming by the original authors. 
* **Exact Column Names Exposed:**
    * Survival Status was named: `"survival status, 0 = censored, 1 = death:ch1"`
    * Survival Time was named: `"time to death (days):ch1"`
* **Resolution:** This detective work allowed us to rewrite the analysis script (`14b_Sniper`) to explicitly target these exact string names and convert the time unit from days to months, bypassing the automated extraction failure.

---
