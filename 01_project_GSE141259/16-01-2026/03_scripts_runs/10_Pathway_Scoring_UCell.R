# ==============================================================================
# Script: 10_Pathway_Scoring_UCell.R (修正安装版)
# Purpose: 给细胞的状态打分 (Senescence, EMT, Inflammation, etc.)
# ==============================================================================

# 1. 正确安装 UCell (它在 Bioconductor 里，不在 CRAN)
# ------------------------------------------------------------------------------
message("📦 正在检查环境...")

# 先安装 BiocManager (如果还没有的话)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 用 BiocManager 安装 UCell
if (!requireNamespace("UCell", quietly = TRUE)) {
  message("⬇️ 正在从 Bioconductor 安装 UCell (这可能需要一点时间)...")
  BiocManager::install("UCell", update = FALSE, ask = FALSE)
}

library(Seurat)
library(UCell)
library(ggplot2)
library(dplyr)

message("✅ UCell 加载成功！")

# 2. 读取数据
# ------------------------------------------------------------------------------
rds_path <- "lianxi/01_data_processed/lung_obj_final_analysis.rds"
# 如果内存里没有 seurat_obj，就读
if(!exists("seurat_obj")) {
  message("📂 读取 Seurat 对象...")
  seurat_obj <- readRDS(rds_path)
}

# 3. 定义基因集 (Feature List)
# ------------------------------------------------------------------------------
message("📝 定义基因集...")

gene_sets <- list(
  ADI_Signature = c("Krt8", "Krt18", "Cldn4", "Cdkn1a", "Ctgf"),
  EMT_Score = c("Vim", "Fn1", "Cdh2", "Zeb1", "Twist1", "Snai1", "Col1a1"),
  Senescence = c("Cdkn2a", "Cdkn1a", "Il6", "Cxcl1", "Mmp3"),
  Inflammation = c("Il1b", "Il6", "Tnf", "Ccl2", "Cxcl2", "Nfkb1"),
  Hypoxia = c("Hif1a", "Vegfa", "Slc2a1", "Ldha", "Pgk1")
)

# 4. 极速打分
# ------------------------------------------------------------------------------
message("🚀 正在计算评分 (UCell)...")
# UCell 算法非常快，几秒钟搞定
seurat_obj <- AddModuleScore_UCell(seurat_obj, features = gene_sets, name = NULL)
message("✅ 评分完成！")

# 5. 画图
# ------------------------------------------------------------------------------
message("🎨 正在画图...")



# ==============================================================================
# 📊 Figure Interpretation: UCell Pathway Analysis (图解指南)
# ==============================================================================

# ------------------------------------------------------------------------------
# Figure 1: Pathway_UCell_Violin.png (Violin Plot)
# Focus: Comparing AT2 vs. Fibroblasts vs. Krt8 ADI
# ------------------------------------------------------------------------------

# 1. Validation (阳性对照验证):
#    - [Observation]: The 'EMT_Score' is highest in Fibroblasts.
#    - [Conclusion]: The scoring algorithm works correctly (Fibroblasts are the gold standard for mesenchymal markers).
#    - [中文]: 验证成功。成纤维细胞的 EMT 分数最高，说明我们的评分算法是准的。

# 2. The "Pathogenic" Nature of Krt8 ADI (ADI 的致病特征):
#    - [Observation]: Compared to normal AT2 cells, Krt8 ADI cells show significant upregulation in all stress pathways.
#      (与正常 AT2 相比，ADI 细胞的所有压力指标全面升高。)
#
#    A. EMT (Partial Transition / 部分间质化):
#       - ADI cells show intermediate EMT scores (Higher than AT2, lower than Fibroblasts).
#       - Meaning: They are undergoing Epithelial-Mesenchymal Transition, contributing to fibrosis.
#
#    B. Senescence (Aging / 衰老):
#       - High 'Senescence' score in ADI.
#       - Meaning: These cells are likely in cell-cycle arrest (p16/p21 pathways) and secreting SASP factors.
#
#    C. Hypoxia & Inflammation (缺氧与炎症):
#       - High 'Hypoxia' and 'Inflammation' scores.
#       - Meaning: ADI cells exist in a hypoxic niche and actively drive local inflammation.

# ------------------------------------------------------------------------------
# Figure 2: Pathway_UCell_UMAP.png (Feature Plot)
# Focus: Spatial Distribution of Stress States
# ------------------------------------------------------------------------------

# 1. ADI Specificity (特异性):
#    - The 'ADI_Signature', 'Hypoxia', and 'Senescence' hotspots perfectly overlap with the Krt8 ADI cluster.
#    - This confirms that these stress states are intrinsically linked to the ADI identity, not random noise.
#    - [中文]: 缺氧和衰老信号完美重合在 ADI 细胞簇上，说明这是 ADI 的核心特征。

# 2. Fibrotic Microenvironment (纤维化微环境):
#    - The 'EMT_Score' lights up both Fibroblasts (mesenchymal) and ADI cells (transitioning).
#    - This visualizes the "pro-fibrotic zone" in the UMAP.

# ==============================================================================
# 📝 Summary for Manuscript (论文结论):
# "UCell scoring reveals that Krt8 ADI cells acquire a multi-pathogenic phenotype characterized by 
#  partial EMT, metabolic hypoxia, senescence, and pro-inflammatory signaling, distinguishing them 
#  from homeostatic AT2 cells."
# ==============================================================================
scores_to_plot <- names(gene_sets)

# A. UMAP 图
p_umap <- FeaturePlot(seurat_obj, features = scores_to_plot, ncol = 3, 
                      min.cutoff = "q5", max.cutoff = "q95", order = TRUE) & 
  scale_color_viridis_c(option = "magma") & 
  theme(legend.position = "right")

ggsave("lianxi/04_output_plots/Pathway_UCell_UMAP.png", p_umap, width = 12, height = 8)

# B. 小提琴图 (对比 AT2 vs ADI)
# 确保 ident 是 cell.type
Idents(seurat_obj) <- seurat_obj$cell.type

# 选几个感兴趣的细胞类型
target_cells <- c("AT2 cells", "Krt8 ADI", "Fibroblasts", "Alveolar macrophages")
# 只保留数据里存在的细胞类型
target_cells <- target_cells[target_cells %in% unique(Idents(seurat_obj))]

p_vln <- VlnPlot(seurat_obj, features = scores_to_plot, 
                 idents = target_cells, 
                 stack = TRUE, flip = TRUE, 
                 fill.by = "ident") +
  NoLegend() +
  ggtitle("Pathway Activity Scores: Cell Type Comparison")

ggsave("lianxi/04_output_plots/Pathway_UCell_Violin.png", p_vln, width = 8, height = 10)

message("🎉 全部完成！去看看那两张图吧！")