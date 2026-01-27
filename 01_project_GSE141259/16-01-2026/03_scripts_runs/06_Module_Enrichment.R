# ==============================================================================
# Script: 06_Module_Enrichment.R
# Purpose: 对 Krt8 ADI 模块基因进行功能富集 (Translation)
# ==============================================================================

library(clusterProfiler)
library(org.Mm.eg.db) # 小鼠数据库
library(dplyr)
library(ggplot2)

# 1. 读取昨天的“嫌疑人名单”
#    (确保这个文件在你的文件夹里)
csv_path <- "lianxi/04_output_plots/Krt8_ADI_Module_Genes.csv"
message("📂 正在读取文件: ", csv_path)

gene_list_df <- read.csv(csv_path)

# 【保险操作】提取基因名
# 因为 write.csv 可能会多出一列序号，我们通常取第二列，或者叫 "x" 的那一列
if("x" %in% colnames(gene_list_df)) {
  target_genes <- gene_list_df$x
} else {
  # 如果没有表头叫 x，就盲猜第二列 (第一列通常是序号 1,2,3...)
  target_genes <- gene_list_df[, 2]
}

message("🕵️‍♀️ 读入基因数量: ", length(target_genes), " 个")
# 看看前几个基因对不对
print(head(target_genes))

# 2. 基因名“翻译” (Symbol -> Entrez ID)
message("🔄 正在转换基因 ID...")
gene_convert <- bitr(target_genes, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db)

# 3. 开始富集分析 (GO Enrichment)
#    审问它们：你们聚在一起搞什么生物学过程 (BP)？
message("🚀 开始运行 GO 富集分析 (可能需要 1-2 分钟)...")
ego <- enrichGO(gene = gene_convert$ENTREZID,
                OrgDb = org.Mm.eg.db,
                ont = "BP",           # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)      # 结果显示基因名

# 4. 画气泡图 (Dotplot)
message("🎨 正在绘图...")
p_dot <- dotplot(ego, showCategory = 15) + 
  ggtitle("Functions of Krt8 ADI Module") +
  theme(axis.text.y = element_text(size = 10)) 

print(p_dot)

# 5. 保存结果
ggsave("lianxi/04_output_plots/Krt8_Module_GO_Enrichment.png", p_dot, width = 8, height = 10)
write.csv(as.data.frame(ego), "lianxi/04_output_plots/Krt8_Module_GO_Table.csv")

message("✅ 审问完成！快看看图里的英文单词是啥？")