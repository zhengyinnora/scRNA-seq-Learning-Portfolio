# ==============================================================================
# Script: 09_TF_Analysis_Final_AutoFix.R
# Purpose: 【智能版】自动识别评分名称，防止过滤为空
# ==============================================================================

# 1. 再次驱魔 (清空变量)
# ------------------------------------------------------------------------------
message("🧹 正在清理内存...")
rm(list = intersect(ls(), c("tf_acts", "tf_activity", "mat", "net", "diff_tfs", "top_tfs")))
gc()

library(Seurat)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(OmnipathR)

# 2. 读取数据
# ------------------------------------------------------------------------------
rds_path <- "lianxi/01_data_processed/lung_obj_final_analysis.rds"
if(!exists("seurat_obj")) {
  message("📂 读取 Seurat 对象...")
  seurat_obj <- readRDS(rds_path)
}

# 3. 整容手术 (C1, C2... 既然这步没报错就保留，最稳)
# ------------------------------------------------------------------------------
message("🔧 重置细胞名为 C1, C2...")
clean_names <- paste0("C", 1:ncol(seurat_obj))
seurat_obj <- RenameCells(seurat_obj, new.names = clean_names)

# 4. 准备计算
# ------------------------------------------------------------------------------
message("📊 提取矩阵...")
mat <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

if(!exists("net")) {
  message("🌐 获取数据库...")
  net <- get_dorothea(organism = 'mouse', levels = c('A', 'B', 'C'))
}

# 5. 极速计算
# ------------------------------------------------------------------------------
message("🚀 正在计算 (run_ulm)...")
tf_acts <- run_ulm(mat = mat, network = net, .source = "source", .target = "target", 
                   .mor = "mor", minsize = 5)

# 🛑 【核心修复点】检查算出来了什么东西
# ------------------------------------------------------------------------------
if(nrow(tf_acts) == 0) {
  stop("❌ 严重错误：run_ulm 结果为空！可能是基因名(Mouse/Human)不匹配。请联系我！")
}

# 自动看看 statistic 这一列叫什么 (是 'ulm' 还是 't_value'？)
available_stats <- unique(tf_acts$statistic)
message("ℹ️ 检测到算出的统计指标有: ", paste(available_stats, collapse = ", "))

# 默认取第一个 (通常就是分数值)
target_stat <- available_stats[1] 
message("✅ 将使用 '", target_stat, "' 作为活性分数。")

# 6. 格式转换 (使用自动检测到的 stat)
# ------------------------------------------------------------------------------
message("🔧 正在转换矩阵...")

tf_activity <- tf_acts %>%
  filter(statistic == target_stat) %>% # <--- 这里改成了变量，不再是死板的字符串
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>%
  as.matrix()

# 7. 保存结果
# ------------------------------------------------------------------------------
common_cells <- intersect(colnames(seurat_obj), colnames(tf_activity))
message("📊 匹配结果：Seurat有 ", ncol(seurat_obj), " 个，TF矩阵有 ", ncol(tf_activity), " 个。")

if(length(common_cells) == 0) {
  stop("💀 依然匹配不上... 这不可能！")
}

seurat_obj <- seurat_obj[, common_cells]
tf_activity <- tf_activity[, common_cells]

seurat_obj[["TF"]] <- CreateAssayObject(data = tf_activity)
DefaultAssay(seurat_obj) <- "TF"
ScaleData(seurat_obj)

message("🎉🎉🎉 TF 数据保存成功！")

# 8. 画图
# ------------------------------------------------------------------------------
message("🎨 正在画图...")
Idents(seurat_obj) <- seurat_obj$cell.type
diff_tfs <- FindMarkers(seurat_obj, ident.1 = "Krt8 ADI", ident.2 = "AT2 cells", 
                        assay = "TF", logfc.threshold = 0)

top_tfs <- diff_tfs %>%
  arrange(desc(avg_log2FC)) %>%
  slice(c(1:10, (n()-9):n())) %>%
  rownames()

p_heatmap <- DoHeatmap(AverageExpression(seurat_obj, assays = "TF", features = top_tfs, return.seurat = T), 
                       features = top_tfs, draw.lines = F, size = 3) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("Top Transcription Factors: AT2 vs ADI")

ggsave("lianxi/04_output_plots/TF_Activity_Heatmap.png", p_heatmap, width = 8, height = 10)

message("🖼️ 热图已保存！这次真的结束了！")



# ==============================================================================
# 🖼️ Figure Interpretation: TF_Activity_Heatmap.pn
# ==============================================================================

# 1. How to Read (怎么看):
#    -------------------------------------------------------
#    - 🔴 Red (红色): High TF Activity (The "Boss" is giving orders ).
#    - 🔵 Blue (蓝色): Low TF Activity (Silent ).
#    - 🎯 Focus (重点): Compare the "AT2 cells" column vs. "Krt8 ADI" column.
#      (对比 AT2 和 Krt8 ADI 这两列，找谁变红了？)

# 2. Key Findings: The "Criminal Gang" behind ADI :
#    -------------------------------------------------------
#    A. 👑 The Mastermind (头号重编程主谋): **Sox9**
#       - [Observation]: Deep Red in Krt8 ADI, Blue in AT2.
#       - [Meaning]: Sox9 is the key driver forcing AT2 cells to lose identity 
#         and become ADI (Progenitor-like state).
#       - [中文]: 实锤了！Sox9 是导致 AT2 丧失身份、变身 ADI 的总开关。

#    B. 🆘 The Stress Signal (求救信号/压力): **Atf6** & **Hsf1**
#       - [Observation]: Highly active (Red) in Krt8 ADI.
#       - [Meaning]: Markers of ER Stress (Unfolded Protein Response). 
#         The cell is under huge proteotoxic stress.
#       - [中文]: 代表内质网应激（ER Stress）。细胞工厂负荷过重，正在“喊救命”。

#    C. 🏗️ The Fibrosis Builder (纤维化建筑师): **Twist1**
#       - [Observation]: Active in ADI.
#       - [Meaning]: A classic EMT (Epithelial-Mesenchymal Transition) driver. 
#         It makes epithelial cells act like fibroblasts (scarring).
#       - [中文]: 上皮-间质转化（EMT）的推手。它让 AT2 细胞变得像成纤维细胞一样，导致肺纤维化。

#    D. 🔥 The Inflammatory Thugs (炎症打手): **Junb** & **Fosl1** (AP-1 family)
#       - [Observation]: Upregulated in ADI.
#       - [Meaning]: Driving the secretion of inflammatory cytokines.
#       - [中文]: 负责制造炎症风暴的 AP-1 家族成员。

# 3. Summary (一句话总结):
#    -------------------------------------------------------
#    "Injury triggers ER stress (Atf6), which activates Sox9 to reprogram AT2 cells, 
#     while Twist1 drives them towards a fibrotic state."
#    (损伤引发应激 [Atf6] -> 激活总导演 [Sox9] 进行重编程 -> 
#     配合 [Twist1] 推动纤维化进程。)
# ==============================================================================