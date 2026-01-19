# ==============================================================================
# Script: 04_Enrichment_analysis.R
# Purpose: Functional enrichment (GO) for Krt8 ADI up-regulated genes.
# ==============================================================================

# 1. 检查并安装包 (如果没有的话)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

# 注意：这是老鼠的数据，必须用 org.Mm.eg.db
if (!require("org.Mm.eg.db", quietly = TRUE))
  BiocManager::install("org.Mm.eg.db")

# 2. 加载包
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

message("包加载完成！准备开始！🚀")

# 1. 安装 clusterProfiler (强制不问问题)
BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)

# 2. 安装 org.Mm.eg.db (强制不问问题)
BiocManager::install("org.Mm.eg.db", update = FALSE, ask = FALSE)


# ==============================================================================
# Step 2: Load Data & Filter (Fixed Version)
# ==============================================================================

# 1. 关键：先把工具箱拿出来！
library(dplyr)             # 提供 %>% 这个符号
library(clusterProfiler)   # 提供 bitr 转换功能
library(org.Mm.eg.db)      # 提供老鼠基因字典

# 2. 读取上次存的 DEG 表格
deg_file <- "lianxi/05_results_tables/DEG_Krt8_vs_AT2.csv"
deg_table <- read.csv(deg_file, row.names = 1)

# 3. 挑选“嫌疑人” (这回肯定能跑通了)
up_genes <- deg_table %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>%
  rownames() 

message("一共找到了 ", length(up_genes), " 个显著上调的基因！")

# 4. 翻译名字 (Symbol -> Entrez ID)
gene_convert <- bitr(up_genes, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db) 

message("成功翻译了 ", nrow(gene_convert), " 个基因 ID。准备进行下一步！")

# ==============================================================================
# Step 3: Run GO Enrichment (Biological Process)
# ==============================================================================

message("开始计算 GO 富集... (可能需要几分钟 ☕️)")

ego <- enrichGO(gene          = gene_convert$ENTREZID, # 输入刚才翻译好的数字ID
                OrgDb         = org.Mm.eg.db,          # 查老鼠库
                ont           = "BP",                  # BP = 生物学过程 (最常用)
                pAdjustMethod = "BH",                  # 统计学矫正方法
                pvalueCutoff  = 0.05,                  # 只要显著的
                qvalueCutoff  = 0.2,
                readable      = TRUE)                  # 结果里把数字ID再翻回英文名，方便你看

message("计算完成！")

# ==============================================================================
# Step 4: Visualization (Dotplot)
# ==============================================================================

# 1. 简单的气泡图 (展示前 15 个最显著的通路)
p_dot <- dotplot(ego, showCategory = 15) + 
  ggtitle("Krt8 ADI Upregulated Pathways (GO:BP)") +
  theme(axis.text.y = element_text(size = 10)) # 字体调大点方便看

print(p_dot)

# 2. 保存图片
ggsave("lianxi/04_output_plots/GO_Enrichment_Dotplot.png", p_dot, width = 8, height = 8)

message("🎉 搞定！图片已保存！快去看看气泡图！")

# ==============================================================================
# [结果解读] GO 富集分析：透视 Krt8 ADI 的“幕后工厂”
# [Result Interpretation] GO Enrichment: Inside the Krt8 ADI Factory
# ------------------------------------------------------------------------------
# 
# 1. 核心发现：疯狂的蛋白质合成 (Hyper-active Protein Synthesis)
#    [Observation]: 
#    - 最显著的通路 (Top terms) 集中在 "cytoplasmic translation" (细胞质翻译) 
#      和 "ribosome biogenesis" (核糖体生成)。
#    [Insight]: 
#    - CN: 这表明 Krt8 ADI 细胞正处于“全功率运转”状态。它们必须疯狂制造大量的蛋白质
#      （如我们之前发现的应激蛋白 S100a6、骨架蛋白 Krt8、伴侣蛋白 Clu）来应对损伤。
#    - EN: The cells are in a hyper-biosynthetic state. They are actively manufacturing 
#      proteins (stress factors, cytoskeletal elements) to survive the injury.
#
# 2. 能量代价：高代谢压力 (High Metabolic Demand)
#    [Observation]:
#    - "aerobic respiration" (有氧呼吸) 和 "oxidative phosphorylation" (氧化磷酸化) 显著富集。
#    [Insight]:
#    - CN: “造蛋白”极其耗能。这群细胞就像一台过载的发动机，正在疯狂燃烧能量 (ATP) 
#      来维持生存。这也暗示了它们处于极高的代谢压力之下。
#    - EN: Protein synthesis is energetically expensive. These terms indicate a surge 
#      in ATP production (oxidative phosphorylation) to fuel the stress response.
#
# 3. 总结 (Conclusion)
#    - CN: Krt8 ADI 细胞虽然停止了分裂（Cdkn1a 高），但绝不是在“休息”。
#      相反，它们内部正在进行一场激烈的“求生战”，透支能量来生产救命物资。
#    - EN: Although these cells are senescent (cycle arrested), they are metabolically 
#      hyper-active, struggling to maintain homeostasis amidst tissue injury.
# ==============================================================================