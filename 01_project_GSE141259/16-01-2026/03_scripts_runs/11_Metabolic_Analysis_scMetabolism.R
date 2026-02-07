# ==============================================================================
# Script: 11_Metabolic_Analysis_UCell_Manual.R
# Purpose: 【绝不报错版】跳过安装，直接用 UCell 算代谢评分
# ==============================================================================

library(Seurat)
library(UCell)
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. 读取数据
# ------------------------------------------------------------------------------
rds_path <- "lianxi/01_data_processed/lung_obj_final_analysis.rds"
if(!exists("seurat_obj")) seurat_obj <- readRDS(rds_path)

# 2. 手动定义核心代谢通路的基因 (Mouse Gene Symbols)
#    这些是 KEGG 数据库里最核心的代谢基因
# ------------------------------------------------------------------------------
message("📝 定义代谢基因集 (直接用硬核列表)...")

metabolic_genes <- list(
  # 🔥 糖酵解 (Glycolysis) - 缺氧/ADI 细胞喜欢用
  Glycolysis = c("Pgk1", "Eno1", "Aldoa", "Ldha", "Gapdh", "Pkm", "Hk2", "Pfkp", "Tpi1", "Gpi1"),
  
  # ⚡ 氧化磷酸化 (OxPhos) - 能量工厂
  OxPhos = c("Cox4i1", "Atp5a1", "Uqcrc2", "Ndufs1", "Ndufb8", "Sdha", "Atp5mc1", "Cycs"),
  
  # 🥥 脂肪酸代谢 (Fatty Acid Metabolism) - 正常 AT2 的看家本领
  Fatty_Acid_Metabolism = c("Fasn", "Acaca", "Scd1", "Acly", "Cpt1a", "Acadm", "Acadl", "Hadh"),
  
  # 🍋 三羧酸循环 (TCA Cycle)
  TCA_Cycle = c("Cs", "Idh1", "Idh2", "Aco2", "Suclg1", "Fh1", "Mdh1", "Mdh2"),
  
  # 🥜 谷胱甘肽代谢 (Glutathione) - 抗氧化/抗压力
  Glutathione = c("Gpx1", "Gpx4", "Gss", "Gsr", "Gstm1", "Gstp1")
)

# 3. 开始打分 (UCell)
# ------------------------------------------------------------------------------
message("🚀 正在计算代谢评分 (UCell)...")
# 不需要 scMetabolism，UCell 算得更准
seurat_obj <- AddModuleScore_UCell(seurat_obj, features = metabolic_genes, name = NULL)

message("✅ 代谢评分完成！")

# 4. 可视化 (气泡图 DotPlot)
# ------------------------------------------------------------------------------
message("🎨 正在绘制代谢气泡图...")

pathways <- names(metabolic_genes)

# DotPlot
p_dot <- DotPlot(seurat_obj, features = pathways, group.by = "cell.type") + 
  coord_flip() + # 翻转
  scale_color_gradient(low = "blue", high = "red") +
  ggtitle("Metabolic Pathway Activity (UCell Score)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("lianxi/04_output_plots/Metabolism_DotPlot_Manual.png", p_dot, width = 9, height = 7)

# 5. 可视化 (FeaturePlot - 重点看糖酵解和脂肪酸)
# ------------------------------------------------------------------------------
message("🎨 正在绘制 UMAP...")

# 看看谁在烧糖 (Glycolysis) vs 谁在烧油 (Fatty_Acid)
p_umap <- FeaturePlot(seurat_obj, features = c("Glycolysis", "Fatty_Acid_Metabolism"), 
                      ncol = 2, min.cutoff = "q10", max.cutoff = "q95", order = T) &
  scale_color_viridis_c(option = "inferno")

ggsave("lianxi/04_output_plots/Metabolism_UMAP_Comparison.png", p_umap, width = 12, height = 6)

# 6. 小提琴图对比 (AT2 vs ADI)
# ------------------------------------------------------------------------------
target_cells <- c("AT2 cells", "Krt8 ADI")
# 检查细胞是否存在
real_targets <- target_cells[target_cells %in% unique(Idents(seurat_obj))]

if(length(real_targets) > 0) {
  p_vln <- VlnPlot(seurat_obj, features = c("Glycolysis", "Fatty_Acid_Metabolism"), 
                   idents = real_targets, pt.size = 0.1)
  ggsave("lianxi/04_output_plots/Metabolism_Violin_AT2vsADI.png", p_vln, width = 8, height = 6)
}

message("🎉 分析完成！")
message("   去看看 Metabolism_DotPlot_Manual.png")
message("   如果 ADI 的 Glycolysis 很高，说明它们正在‘癌变’式生存！")



# ==============================================================================
# 📊 Figure Interpretation: Metabolic Reprogramming Analysis (图解指南)
# ==============================================================================

# ------------------------------------------------------------------------------
# Figure 1: Metabolism_DotPlot_Manual.png (DotPlot)
# Focus: The "Metabolic Switch" between AT2 and ADI
# ------------------------------------------------------------------------------

# 1. AT2 Identity = Lipid Metabolism (AT2 的看家本领: 脂代谢):
#    - [Observation]: The 'Fatty_Acid_Metabolism' dot is large and deep red in 'AT2 cells' and 'Activated AT2'.
#    - [Biological Meaning]: AT2 cells require high lipid metabolism to synthesize pulmonary surfactant (dipalmitoylphosphatidylcholine). 
#      This confirms their functional homeostatic state.
#    - [中文]: AT2 细胞高表达脂肪酸代谢通路，这是它们合成肺表面活性物质（脂质）的基础，代表正常功能状态。

# 2. ADI Reprogramming = Glycolysis (ADI 的重编程: 糖酵解):
#    - [Observation]: In 'Krt8 ADI' cells, 'Fatty_Acid_Metabolism' drops (turns blue/small), 
#      while 'Glycolysis' becomes highly active (Red).
#    - [Biological Meaning]: This indicates a "Warburg-like" metabolic shift. 
#      ADI cells abandon surfactant production and switch to glycolysis to survive in the hypoxic injury niche 
#      and fuel rapid plasticity/proliferation.
#    - [中文]: 代谢发生大逆转。ADI 细胞放弃了耗氧的脂代谢，转而依赖糖酵解。这是为了适应缺氧环境（Hypoxia），并为细胞变形提供快速能量。

# 3. Oxidative Stress Response (氧化应激):
#    - [Observation]: 'Glutathione' metabolism is upregulated in ADI.
#    - [Meaning]: Consistent with the UCell 'Senescence' score, ADI cells are actively fighting ROS (Reactive Oxygen Species) stress.

# ------------------------------------------------------------------------------
# Figure 2: Metabolism_UMAP_Comparison.png (FeaturePlot)
# Focus: Spatial Exclusion of Metabolic States
# ------------------------------------------------------------------------------

# 1. Mutually Exclusive Domains (互斥的代谢区域):
#    - The UMAP clearly shows that cells high in 'Glycolysis' (left plot, bright spots) 
#      do NOT overlap with cells high in 'Fatty_Acid_Metabolism' (right plot).
#    - [Conclusion]: The loss of lipid metabolism is a hallmark of the AT2-to-ADI transition.
#      (脂代谢的丢失是 AT2 向 ADI 转变的关键标志。)

# ==============================================================================
# 📝 Summary for Manuscript (论文结论):
# "Metabolic scoring reveals a profound reprogramming during the AT2-to-ADI transition: 
#  cells downregulate surfactant-associated fatty acid metabolism and upregulate glycolysis, 
#  likely as an adaptive strategy to the hypoxic injury microenvironment."
# ==============================================================================