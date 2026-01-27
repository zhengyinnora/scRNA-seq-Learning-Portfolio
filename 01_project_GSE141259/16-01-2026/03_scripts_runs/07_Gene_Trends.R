# ==============================================================================
# Script: 07_Gene_Trends.R
# Purpose: 画出明星基因随拟时序变化的趋势线 (Line Plots)
# ==============================================================================

library(monocle3)
library(ggplot2)
library(dplyr)

# 1. 既然是新脚本，必须重新加载数据
#    (读取昨天存好的那个终极 CDS 对象)
cds_path <- "lianxi/01_data_processed/lung_monocle_final.rds"
if(file.exists(cds_path)) {
  message("📂 正在读取 CDS 对象...")
  cds <- readRDS(cds_path)
} else {
  stop("❌ 找不到 lung_monocle_final.rds！请检查路径或先运行 05 脚本保存数据！")
}

# 2. 定义你要看的“明星基因”
#    我们选几个经典的：
#    - Sftpc, Sftpb (AT2 的看家基因，理论上应该下降)
#    - Krt8, Krt19, Lgals3 (你抓到的坏蛋，理论上应该上升)
#    - Mki67 (看看有没有细胞在偷偷分裂？)
target_genes <- c("Sftpc", "Sftpb", "Krt8", "Krt19", "Lgals3", "Mki67")

# 检查一下这些基因在不在数据里
valid_genes <- target_genes[target_genes %in% rowData(cds)$gene_short_name]

# 3. 画趋势图 (Feature Plot along Pseudotime)
message("🎨 正在绘制基因趋势图...")

# plot_genes_in_pseudotime 是 Monocle3 专门画这种图的函数
# min_expr: 过滤掉表达量太低的噪音
p_trends <- plot_genes_in_pseudotime(cds[valid_genes,],
                                     color_cells_by = "cell_type", # 点的颜色按细胞类型分
                                     min_expr = 0.5, 
                                     ncol = 2) +                   # 排成 2 列
  scale_color_manual(values = c("AT2 cells" = "#440154FF",     # 紫色
                                "Krt8 ADI" = "#FDE725FF")) +   # 黄色
  ggtitle("Gene Expression Trends along Pseudotime")

print(p_trends)

# 4. 保存
ggsave("lianxi/04_output_plots/Monocle3_Gene_Trends_LinePlot.png", p_trends, width = 8, height = 10)

message("✅ 趋势图画好了！快看看 Sftpc 是不是在跌，Krt8 是不是在涨？")


# ==============================================================================
# [Result Interpretation: Gene Expression Dynamics (Refined)]
# ==============================================================================
# 1. Loss of Identity (Sftpc/Sftpb):
#    - Show a consistent and sharp decline along pseudotime, confirming the
#      gradual loss of AT2 lineage identity.
#
# 2. Heterogeneous Gain of ADI Markers (Krt8/Krt19):
#    - Although the smoothed trend line remains low due to sparsity (dropouts),
#      a distinct population of cells at the terminus (Yellow/ADI phase) 
#      shows high expression levels (dots > 3.0), contrasting with the
#      complete absence in the early phase.
#
# 3. Transient Activation (Lgals3):
#    - Lgals3 exhibits the most robust dynamic pattern, showing a clear 
#      "bell-shaped" upregulation during the intermediate transition phase,
#      identifying it as a high-confidence driver of the transdifferentiation.
# ==============================================================================