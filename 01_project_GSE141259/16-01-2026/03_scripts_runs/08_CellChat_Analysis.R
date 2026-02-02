# ==============================================================================
# Script: 08_CellChat_Analysis.R
# Purpose: 安装 CellChat 并准备环境
# ==============================================================================
# step1
# 1. 安装 devtools (如果还没有的话)
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# 2. 安装 NMF (CellChat 依赖这个做模式识别)
#    这个包有时候很难装，如果报错，通常是因为缺系统库，先试着装
if (!requireNamespace("NMF", quietly = TRUE))
  install.packages("NMF")

# 3. 安装 Circlize and ComplexHeatmap (画图用的)
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}

# 4. 安装 CellChat (主角)
#    直接从 GitHub 装最新版
if (!requireNamespace("CellChat", quietly = TRUE)) {
  message("🚀 正在下载并安装 CellChat，这可能需要一点时间...")
  devtools::install_github("sqjin/CellChat")
}

library(CellChat)
library(ggplot2)
library(dplyr)
library(Seurat)

message("✅ CellChat 环境加载成功！第一关通关！")


# ==============================================================================
# step2：Data Prep: Seurat -> CellChat
# ==============================================================================

# 1. 读取之前的 Seurat 对象
#    注意：我们要用包含所有细胞的对象，不仅仅是 AT2/ADI，
#    因为细胞通讯看的是“别人”怎么跟 AT2 说话。
#    如果你只有 lung_obj_final_analysis.rds (子集)，也没关系，先用这个练手。
#    如果有更大的包含免疫细胞的对象最好，没有就先跑当前的。
rds_path <- "lianxi/01_data_processed/lung_obj_final_analysis.rds" 

if(file.exists(rds_path)) {
  seurat_obj <- readRDS(rds_path)
  message("📂 数据加载成功！")
} else {
  stop("❌ 找不到 Seurat 对象文件！")
}

# 2. 提取数据矩阵 (Data Matrix)
#    CellChat 需要标准化的数据 (NormalizeData 之后的数据)
data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

# 3. 提取细胞类型 (Meta Data)
#    确保你的 Seurat 对象里有 'cell.type' 这一列
#    如果有其他名字（比如 cell_type），请在这里修改
meta <- seurat_obj@meta.data
# 检查一下列名，防止报错
if("cell.type" %in% colnames(meta)) {
  cell_labels <- meta$cell.type
} else {
  # 假如列名叫其他的，比如 active.ident
  cell_labels <- Idents(seurat_obj)
}

# 把细胞类型变成一个数据框
meta_df <- data.frame(labels = cell_labels, row.names = rownames(meta))

# 4. 创建 CellChat 对象
message("🔨 正在创建 CellChat 对象...")
cellchat <- createCellChat(object = data.input, meta = meta_df, group.by = "labels")

# 5. 设置数据库 (Human or Mouse)
#    你的数据是小鼠，所以用 CellChatDB.mouse
message("📚 加载受体-配体数据库 (Mouse)...")
CellChatDB <- CellChatDB.mouse 
# 使用全部数据库 (包含 Secreted Signaling, ECM-Receptor, Cell-Cell Contact)
cellchat@DB <- CellChatDB

# 6. 预处理 (Subset data & Identify Over-expressed genes)
#    这一步是把表达量算一下，找找哪些受体配体是高表达的
message("⚙️ 正在进行预处理计算 (Pre-processing)...")
cellchat <- subsetData(cellchat) # 这一步是必要的
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

message("✅ 数据准备完成！可以开始推断通讯网络了！")


# ==============================================================================
# step3： Inference: 推断细胞通讯网络
# ==============================================================================

message("🧠 开始推断通讯概率 (Compute Communication Probability)...")

# 1. 计算通讯概率
#    type = "triMean" 是推荐的计算方法，比较稳健
cellchat <- computeCommunProb(cellchat, type = "triMean")

# 2. 过滤掉细胞太少的通讯 (比如某种细胞只有几个，算出来的结果不可信)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 3. 提取推断出的通讯网络 (这是最后画图用的数据框)
df.net <- subsetCommunication(cellchat)

# 4. 计算通路层面的通讯 (Compute Pathway Level)
#    比如把 "Tgfb1-Tgfbr1" 和 "Tgfb2-Tgfbr2" 都归纳为 "TGFb" 通路
cellchat <- computeCommunProbPathway(cellchat)

# 5. 计算整合网络 (Aggregated Network)
cellchat <- aggregateNet(cellchat)

message("🎉 计算全部完成！CellChat 对象已就绪！")
# 保存一下，防止崩了白跑
saveRDS(cellchat, "lianxi/01_data_processed/lung_cellchat.rds")

# ==============================================================================
# step 4： Visualization: 看看谁在跟谁说话？
# ==============================================================================

# 画一张所有细胞之间互动的数量图
groupSize <- as.numeric(table(cellchat@idents)) # 细胞数量

par(mfrow = c(1,2), xpd=TRUE) # 设置画布

# 图1：互动数量 (Number of interactions)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")

# 图2：互动强度 (Interaction weights/strength)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights")

message("👀 快看 Plots 面板！那个圆圆的图就是细胞社交网络！")



# ==============================================================================
# step 5 ：Save Plots: 保存高清大图 (Base R 专用法)
# ==============================================================================

message("💾 正在保存全景通讯图...")

# 方法 A: 保存为 PDF (推荐！矢量图，无限放大不模糊，文字可以选中)
# ------------------------------------------------------------------
pdf("lianxi/04_output_plots/CellChat_Circle_All.pdf", width = 12, height = 8)

par(mfrow = c(1,2), xpd=TRUE) # 重新设置画布
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights")

dev.off() # 关掉设备 (这一步这才是真正的保存！)


# 方法 B: 保存为 PNG (如果不方便看 PDF，就存这个超大分辨率的图片)
# ------------------------------------------------------------------
# ==============================================================================
# Save Plots Separately: 分别保存，杜绝覆盖
# ==============================================================================

message("💾 正在分别保存高清大图...")

# 1. 保存第一张：数量图 (Number)
# ------------------------------------------------------------------
png("lianxi/04_output_plots/CellChat_Net_Number.png", width = 1500, height = 1500, res = 300)
par(xpd=TRUE) # 允许文字溢出边界，防止被切掉
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
dev.off() # 📸 咔嚓，保存第一张

# 2. 保存第二张：强度图 (Weights)
# ------------------------------------------------------------------
png("lianxi/04_output_plots/CellChat_Net_Weights.png", width = 1500, height = 1500, res = 300)
par(xpd=TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights")
dev.off() # 📸 咔嚓，保存第二张

message("✅ 搞定！现在文件夹里应该有两个文件了，互不干扰！")

# ==============================================================================
# [Result Interpretation: Global Communication Network]
# ==============================================================================
# 1. Overview (总体概览):
#    - These Chord Diagrams visualize the aggregated cell-cell communication network
#      across all cell groups in the dataset.
#    - Nodes (Circle edge): Represent different cell types (e.g., AT2, ADI, Macrophages).
#    - Edges (Lines): Represent inferred communication links based on Ligand-Receptor co-expression.
#
# 2. Plot A: Number of Interactions (数量图):
#    - Visualizes the *Frequency* of communication.
#    - Thicker lines indicate that Cell A and Cell B share a higher number of 
#      active Ligand-Receptor pairs.
#    - "Who talks the most?" (Quantity).
#
# 3. Plot B: Interaction Weights (强度图):
#    - Visualizes the *Strength* of communication (Probability).
#    - Thicker lines indicate high expression levels of the signaling molecules
#      involved in the interaction.
#    - "Who shouts the loudest?" (Quality/Intensity).
#
# 4. Biological Insight (生物学意义):
#    - The dense network connectivity (the "Hairball") suggests a highly coordinated
#      multicellular response to lung injury.
#    - It confirms that the AT2-to-ADI transition does not happen in isolation but
#      is likely regulated by extensive cross-talk with the microenvironment 
#      (Immune cells, Fibroblasts, etc.).
# ==============================================================================


# ==============================================================================
# Part 2: Zoom In - 定向分析 (Bubble Plots)
# ==============================================================================

# 1. 准备工作：确认你的细胞名字
#    (我从你刚才的图里看到了 "AT2 cells" 和 "Krt8 ADI"，所以直接用这俩名字)
#    如果不报错，说明名字是对的。
target_AT2 <- "AT2 cells"   # 受害者
source_ADI <- "Krt8 ADI"    # 嫌疑人

# 2. 场景一：谁在霸凌 AT2？(Incoming Signals to AT2)
#    我们想看所有细胞发给 AT2 的信号
message("🔍 正在分析 AT2 收到的信号 (Incoming)...")

# sources.use = NULL 代表“所有人”
# targets.use = target_AT2 代表“发给 AT2”
# remove.isolate = FALSE 保证即使信号弱也显示出来
p_bubble_incoming <- netVisual_bubble(cellchat, 
                                      sources.use = NULL, 
                                      targets.use = target_AT2, 
                                      remove.isolate = FALSE) +
  ggtitle("Signals Received by AT2 cells") +
  coord_flip() # 翻转一下，方便阅读基因名

# 3. 场景二：Krt8 ADI 在向外散布什么？(Outgoing Signals from ADI)
#    我们想看 ADI 发给别人的信号
message("📢 正在分析 ADI 发出的信号 (Outgoing)...")

# sources.use = source_ADI 代表“来自 ADI”
# targets.use = NULL 代表“发给所有人”
p_bubble_outgoing <- netVisual_bubble(cellchat, 
                                      sources.use = source_ADI, 
                                      targets.use = NULL, 
                                      remove.isolate = FALSE) +
  ggtitle("Signals Sent by Krt8 ADI") +
  coord_flip()

# 4. 保存这两张清爽的图
message("💾 正在保存气泡图...")
ggsave("lianxi/04_output_plots/CellChat_Incoming_AT2.png", 
       p_bubble_incoming, width = 10, height = 12) # 这个图可能比较长，高设大点
ggsave("lianxi/04_output_plots/CellChat_Outgoing_ADI.png", 
       p_bubble_outgoing, width = 10, height = 12)

# 5. 显示出来看看
print(p_bubble_incoming)

# ==============================================================================
# [Biological Interpretation: Ligand-Receptor Specifics]
# ==============================================================================
# 1. Incoming Signals to AT2 (Environment Stress):
#    - Observation: AT2 cells show strong reception of the "Fn1 - Sdc4" signaling pair.
#    - Senders: Primarily Fibroblasts, Myofibroblasts, and Fn1+ Macrophages.
#    - Meaning: The accumulation of Extracellular Matrix (Fibronectin) acts as a 
#      mechanical stress signal received by AT2 cells, likely driving their 
#      loss of identity.
#
# 2. Outgoing Signals from Krt8 ADI (Pro-inflammatory Feedback):
#    - Observation: Krt8 ADI cells strongly upregulate "Mif" and "Spp1" ligands.
#    - Receivers: Macrophages (AM, M2) and Proliferating cells via "Cd74" and "Cd44" receptors.
#    - Meaning: The ADI state actively perpetuates inflammation by recruiting and 
#      activating immune cells (Mif/Spp1 axes), creating a pathogenic loop.
# ==============================================================================