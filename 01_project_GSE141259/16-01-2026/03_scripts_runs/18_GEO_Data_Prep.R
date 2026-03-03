# ==============================================================================
# Script: 18_GEO_Data_Prep.R (修复版 - 探针转换)
# Purpose: 破解芯片探针乱码，提取精准的机器学习矩阵
# ==============================================================================

library(GEOquery)
library(dplyr)

message("🌐 1. 正在加载 GEO 数据 (带着注释文件)...")
# AnnotGPL = TRUE 会顺带把官方的“探针翻译字典”下下来
gse <- getGEO("GSE32537", GSEMatrix = TRUE, AnnotGPL = TRUE)
expr_data <- exprs(gse[[1]])
pheno_data <- pData(gse[[1]])
f_data <- fData(gse[[1]]) # 这里面装着探针和基因的对应关系

message("🧬 2. 正在把讨厌的探针ID翻译成真正的基因名...")
# 提取官方注释里的 Gene symbol 列
gene_symbols <- f_data$`Gene symbol`

# 芯片数据经常一个探针对应多个基因(用 /// 隔开)，我们只取第一个最准确的
gene_symbols <- sapply(strsplit(as.character(gene_symbols), "///"), `[`, 1)
gene_symbols <- trimws(gene_symbols)

# 过滤掉那些没有对应基因名（空白）的垃圾探针
valid_idx <- which(gene_symbols != "" & !is.na(gene_symbols))
expr_data <- expr_data[valid_idx, ]
rownames(expr_data) <- gene_symbols[valid_idx]

# 如果一个基因有多个探针，计算平均值合并（这也是标准操作）
expr_data <- as.data.frame(expr_data) %>%
  mutate(Gene = rownames(.)) %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean)) %>%
  tibble::column_to_rownames("Gene")

message("🧹 3. 正在提取临床分组信息...")
clinical_info <- data.frame(
  Sample = rownames(pheno_data),
  Status = ifelse(grepl("idiopathic pulmonary fibrosis", pheno_data$title, ignore.case = TRUE), "IPF", "Control")
)

message("🎯 4. 正在导入我们的跨物种靶点并提取交集...")
my_targets <- read.csv("lianxi/04_output_plots/17_Drug_Repurposing_Candidates.csv") %>% 
  pull(gene_name) %>% 
  unique()

# 这次是真正的基因名和基因名的碰撞了！
common_genes <- intersect(rownames(expr_data), my_targets)
message(paste("🎉 完美！成功在临床队列中匹配到了", length(common_genes), "个核心致病基因！"))

# 提取子集并转置给机器学习用
ml_data <- as.data.frame(t(expr_data[common_genes, ]))
ml_data$Status <- clinical_info$Status

message("💾 5. 正在保存完美的机器学习输入数据...")
saveRDS(ml_data, "lianxi/01_data_processed/18_ML_Input_Data.rds")
message("✅ 18 步修复通关！")