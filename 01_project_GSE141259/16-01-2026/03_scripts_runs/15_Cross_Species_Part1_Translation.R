# ==============================================================================
# Script: 15_Cross_Species_Part1_Translation_Magic.R
# Purpose: 跨物种基因翻译 (稀疏矩阵乘法终极优化版，0.5秒极速合并)
# ==============================================================================

library(Seurat)
library(dplyr)
library(nichenetr) 
library(Matrix)

# 1. 加载小鼠 Seurat 对象
# ------------------------------------------------------------------------------
message("📂 正在加载小鼠 Seurat 对象...")
rds_path <- "lianxi/01_data_processed/lung_obj_final_analysis.rds"
if(!exists("seurat_obj")) seurat_obj <- readRDS(rds_path)

# 2. 提取并翻译基因
# ------------------------------------------------------------------------------
message("🧬 正在提取矩阵并翻译基因名...")
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
mouse_genes <- rownames(counts)
human_genes <- convert_mouse_to_human_symbols(mouse_genes)

# 过滤掉无法翻译的基因
valid_idx <- !is.na(human_genes)
counts_valid <- counts[valid_idx, ]
human_genes_valid <- human_genes[valid_idx]

# 3. 终极魔法：稀疏矩阵相乘实现极速合并！⚡
# ------------------------------------------------------------------------------
message("⚡ 正在执行稀疏矩阵乘法合并 (见证奇迹的 0.5 秒)...")

# 找出所有独立的人类基因名
unique_human_genes <- unique(human_genes_valid)

# 构建一个"映射矩阵" (只记录 1，保留极其稀疏的状态)
mapping_matrix <- sparseMatrix(
  i = match(human_genes_valid, unique_human_genes), # 行：人类基因索引
  j = seq_along(human_genes_valid),                 # 列：原小鼠基因索引
  x = 1,                                            # 值全为 1
  dims = c(length(unique_human_genes), length(human_genes_valid)),
  dimnames = list(unique_human_genes, rownames(counts_valid))
)

# 魔法发生：(人类基因 x 小鼠基因) %*% (小鼠基因 x 细胞) = (人类基因 x 细胞)
counts_humanized_sparse <- mapping_matrix %*% counts_valid

# 4. 构建带有人类签证的新对象
# ------------------------------------------------------------------------------
message("🏗️ 正在构建全新 Seurat 对象...")
mouse_humanized <- CreateSeuratObject(counts = counts_humanized_sparse, 
                                      meta.data = seurat_obj@meta.data)

mouse_humanized$Species <- "Mouse"
mouse_humanized$Dataset <- "Bleomycin_Model"

message("🔄 正在重新跑基础降维 (稍微需要十几秒)...")
mouse_humanized <- NormalizeData(mouse_humanized) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData()

# 5. 保存
# ------------------------------------------------------------------------------
message("💾 正在保存...")
save_path <- "lianxi/01_data_processed/mouse_humanized_for_integration.rds"
saveRDS(mouse_humanized, save_path)

message("🎉 搞定！不仅风扇没转，而且快到飞起！")
message(paste("✅ 成功保留了", nrow(mouse_humanized), "个人类同源基因。"))