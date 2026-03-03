# ==============================================================================
# Script: 17_Drug_Repurposing.R (自给自足终极版 - 列名修复)
# Purpose: 直接从 RDS 提取靶点并直连数据库找药，已匹配 2026 年最新官网列名
# ==============================================================================

library(Seurat)
library(dplyr)
library(readr)

message("📂 1. 正在加载跨物种 RDS 对象...")
merged_obj <- readRDS("lianxi/04_output_plots/Cross_Species_Integrated_Final.rds")
DefaultAssay(merged_obj) <- "RNA"

message("🔍 2. 正在提取跨物种核心靶点...")
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 200)
target_genes <- VariableFeatures(merged_obj)
target_genes <- toupper(target_genes)

message("✅ 成功提取到 ", length(target_genes), " 个高潜力靶点！")

message("🌐 3. 正在直连 DGIdb 官网下载最新药物数据库...")
dgidb_url <- "https://www.dgidb.org/data/latest/interactions.tsv"
dgidb_data <- read_tsv(dgidb_url, show_col_types = FALSE)

message("💊 4. 正在比对靶点并寻找潜在的抑制剂药物...")
# 核心修复区：使用了你截图里提取出的最新列名
candidate_drugs <- dgidb_data %>%
  filter(gene_name %in% target_genes) %>%
  filter(grepl("inhibitor|antagonist|blocker|suppressor|antibody", interaction_type, ignore.case = TRUE)) %>%
  select(gene_name, drug_name, interaction_type, interaction_source_db_name, interaction_score) %>%
  distinct() %>%
  arrange(desc(interaction_score), gene_name) # 按药物可靠性得分降序排列

write.csv(candidate_drugs, "lianxi/04_output_plots/17_Drug_Repurposing_Candidates.csv", row.names = FALSE)

message("🎉 终极整合版彻底跑通！快去 lianxi/04_output_plots/ 文件夹查看 CSV 吧！")