# ==============================================================================
# Project: GSE141259 Mouse Lung (Strunz et al. 2020)
# Purpose: Load 10x Genomics Data & Basic QC
# Date: 2025-12-12
# ==============================================================================

# 1. 加载包 -------------------------------------------------------------------
library(Seurat)
library(tidyverse)

# 2. 读取数据 (Read 10x Data) ------------------------------------------------
# 使用 ReadMtx 读取三个文件
# 这里的路径必须跟你文件夹里的一模一样！

print("正在读取数据... (可能需要几秒钟)")

raw_counts <- ReadMtx(
  mtx = "01_project_GSE141259/00_data_raw/GSE141259_WholeLung_rawcounts.mtx.gz",
  cells = "01_project_GSE141259/00_data_raw/GSE141259_WholeLung_barcodes.txt.gz",
  features = "01_project_GSE141259/00_data_raw/GSE141259_WholeLung_genes.txt.gz",
  feature.column = 1  # 👈 关键就是加了这一行！
)

print("读取成功！正在创建 Seurat 对象...")

# 3. 创建 Seurat 对象 (Create Object) ----------------------------------------
strunz_obj <- CreateSeuratObject(
  counts = raw_counts,
  project = "Strunz_Lung",
  min.cells = 3,       # 过滤掉极少表达的基因
  min.features = 200   # 过滤掉质量太差的细胞
)

# 4. 看看长啥样 (Inspect) ----------------------------------------------------
print(strunz_obj) 
# 期待输出：你会看到它是多少个细胞 (samples) x 多少个基因 (features)

print("恭喜！数据已成功载入内存！下一步可以做 QC 了！🎉")

# ==============================================================================
# Phase 2: 处理与出图 (Processing & Plotting)
# ==============================================================================

print("🚀 Phase 2 启动！正在创建 Seurat 对象...")

# 1. 创建对象 ------------------------------------------------------------------
strunz_obj <- CreateSeuratObject(
  counts = raw_counts,
  project = "Strunz_Lung",
  min.cells = 3,
  min.features = 200
)

# 🧹 内存清理小技巧：原来的大矩阵没用了，删掉它释放内存
rm(raw_counts)
gc()

print(strunz_obj) # 看看咱们有多少兵马

# 2. 抽样 (Subsampling) - 关键一步！--------------------------------------------
# 真实数据有近 3 万个细胞，跑全套可能要半小时。
# 为了周五的快乐，我们随机抽 5000 个先看看长啥样！
print("🎲 正在随机抽取 5000 个细胞练手...")
set.seed(123)
small_obj <- subset(strunz_obj, cells = sample(Cells(strunz_obj), 5000))

# 3. 标准流程一把梭 (Standard Pipeline) ----------------------------------------
print("⚙️ 流水线开启：标准化 -> 降维 -> 聚类...")
small_obj <- small_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:15) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = 0.5)

# 4. 见证奇迹的时刻 (Plotting) -------------------------------------------------
print("🎨 正在作画...")

# 图A: 分群图 (看看有多少种细胞)
p1 <- DimPlot(small_obj, label = TRUE) + ggtitle("Strunz Lung (5k Subset)")

# 图B: 找找巨噬细胞 (AM Marker: Adgre1)
print("🕵️‍♀️ 嫌疑人已锁定！正在绘制巨噬细胞地图...")

# 改用 "Emr1" (Adgre1 的旧名字)
# 顺便把 Marco 也画出来，双重确认！
p2 <- FeaturePlot(small_obj, features = c("Emr1", "Marco", "Cd68"), ncol = 3)

# 显示图片
p2

# 保存这幅胜利的画作
ggsave("01_project_GSE141259/04_output_plots/Final_AM_Markers_Validated.png", plot = p2, width = 12, height = 4)

print("🎉 搞定！Nora，快看图！那些亮起来的紫色点点，就是你的 AM！")
# 显示图片
p1
p2

# 5. 保存战果 ------------------------------------------------------------------
ggsave("01_project_GSE141259/04_output_plots/Final_UMAP.png", plot = p1)
ggsave("01_project_GSE141259/04_output_plots/Final_AM_Marker.png", plot = p2)

print("🎉 大功告成！Nora，你可以下班了！")

# 1. 看看前 10 个基因名字长啥样？(是 Adgre1 还是 ENSMUSG...?)
head(rownames(small_obj))

# 2. 帮我全网通缉 "Adgre1" (忽略大小写搜索)
# 如果它叫 ADGRE1 或者 Adgre1-a 啥的，这行代码能把它抓出来
grep("Adgre1", rownames(small_obj), value = TRUE, ignore.case = TRUE)

print("正在全网通缉 AM 标记基因...")

# 1. 试试曾用名 Emr1 (Adgre1 的旧名字)
grep("Emr1", rownames(small_obj), value = TRUE, ignore.case = TRUE)

# 2. 试试其他的巨噬细胞铁杆 Marker (总得有一个在吧！)
# Siglecf (肺泡巨噬细胞特异)
grep("Siglecf", rownames(small_obj), value = TRUE, ignore.case = TRUE)

# Marco (也是 AM 的标志)
grep("Marco", rownames(small_obj), value = TRUE, ignore.case = TRUE)

# Cd68 (所有巨噬细胞都有)
grep("Cd68", rownames(small_obj), value = TRUE, ignore.case = TRUE)

# ==============================================================================
# 5. 终极绘图 (Final Plotting)
# ==============================================================================

print("🎨 正在绘制巨噬细胞地图 (Emr1, Marco, Cd68)...")

# 这里的 Emr1 就是 Adgre1 的真身！
# 我们把这三个巨噬细胞的标志物画在一起
p_final <- FeaturePlot(small_obj, features = c("Emr1", "Marco", "Cd68"), ncol = 3)

# 显示图片
p_final

# 保存图片
ggsave("01_project_GSE141259/04_output_plots/Final_AM_Markers_Validated.png", plot = p_final, width = 12, height = 4)

print("🎉 恭喜！Strunz 2020 数据复现成功！快去 output 文件夹看图！")


print("🔎 正在侦察其他细胞类型...")

# 常用的小鼠肺部 Marker：
# Epcam = 上皮细胞 (肺泡上皮)
# Col1a1 = 成纤维细胞 (Fibroblasts, 负责纤维化的坏蛋)
# Cd3e = T 细胞
# Cd79a = B 细胞
# Pecam1 (CD31) = 内皮细胞 (血管)

p_neighbors <- FeaturePlot(small_obj, 
                           features = c("Epcam", "Col1a1", "Cd3e", "Cd79a", "Pecam1"), 
                           ncol = 3)

p_neighbors
ggsave("01_project_GSE141259/04_output_plots/Neighbors.png", plot = p_neighbors, width = 12, height = 8)


print("🎯 正在为 Nora 的课题筛选趋化因子受体...")

# 这些是著名的单核/巨噬细胞迁移受体：
# Ccr2 (招募单核细胞的关键)
# Cx3cr1 (驻留型巨噬细胞)
# Ccr5, Cxcr4 (常见的迁移受体)
# Itgam (CD11b, 也是迁移相关的整合素)

# 我们用气泡图 (DotPlot) 来一次性看清它们在不同群里的表达
p_screen <- DotPlot(small_obj, 
                    features = c("Emr1", "Marco", "Ccr2", "Cx3cr1", "Ccr5", "Cxcr4", "Itgam")) + 
  RotatedAxis() + # 让横坐标文字斜过来，防止重叠
  ggtitle("Potential Migration Targets Screen")

p_screen
ggsave("01_project_GSE141259/04_output_plots/Chemokine_Screen.png", plot = p_screen, width = 8, height = 5)

print("💾 正在存档，请勿关闭 RStudio...")

# 保存处理好的小对象 (5000个细胞版)
saveRDS(small_obj, "01_project_GSE141259/GSE141259_small_processed.rds")

print("✅ 存档完毕！下次直接 readRDS 就能接着玩，不用从头跑了！")



print("🎯 正在响应 Nora 的号召，添加 FPR 家族...")

# 这次我们把名单列全一点！
# 1. 身份卡: Emr1, Marco
# 2. 细菌雷达 (FPRs): Fpr1, Fpr2
# 3. 迁移受体 (CCRs/CXCRs): Ccr2, Ccr5, Cxcr4, Cx3cr1

genes_to_check <- c("Emr1", "Marco", 
                    "Fpr1", "Fpr2",        # 👈 这里！加上了！
                    "Ccr2", "Ccr5", "Cxcr4", "Cx3cr1")

# 画气泡图
p_screen_v2 <- DotPlot(small_obj, features = genes_to_check) + 
  RotatedAxis() +
  ggtitle("Receptors Screen (Including FPRs)") +
  scale_color_gradient(low = "grey90", high = "red") # 换个醒目的红色

p_screen_v2

# 保存
ggsave("01_project_GSE141259/04_output_plots/Chemokine_Screen_with_FPR.png", plot = p_screen_v2, width = 9, height = 5)




print("🏷️ 正在给细胞贴标签 (Annotating)...")

# 1. 定义 6 大门派的“掌门人” (Marker Genes)
# 这一步是为了看清楚哪个 Cluster 是哪种细胞
markers <- c(
  "Emr1",   # AM (巨噬细胞)
  "Epcam",  # Epithelial (上皮: AT1/AT2)
  "Col1a1", # Fibroblast (成纤维)
  "Pecam1", # Endothelial (血管)
  "Cd3e",   # T cells
  "Cd79a",  # B cells
  "S100a8"  # Neutrophils (中性粒细胞)
)

# 2. 先画一张图，看看 0, 1, 2 分别是谁？
# 这张图的横坐标是数字，纵坐标是基因
p_check <- DotPlot(small_obj, features = markers) + 
  RotatedAxis() +
  ggtitle("Who is who?")

p_check
ggsave("01_project_GSE141259/04_output_plots/Identity_Check.png", plot = p_check, width = 8, height = 5)

print("👀 请看刚才生成的 Identity_Check 图！")
print("👉 哪个 Cluster 的 Emr1 最红，那个 Cluster 就是 AM！")
print("👉 哪个 Cluster 的 Epcam 最红，那个 Cluster 就是上皮！")


print("⚖️ 正在进行脱靶效应排查 (Specificity Check)...")

# 组合拳：身份基因 + 你的目标受体
# 前面是身份 (用来定位细胞类型)，后面是靶点 (用来检查干扰)
check_list <- c(
  # --- 身份区 ---
  "Emr1", "Marco",    # AM
  "Col1a1",           # 成纤维 (干扰项1)
  "Epcam",            # 上皮 (干扰项2)
  "Pecam1",           # 血管 (干扰项3)
  "Cd3e",             # T细胞 (干扰项4)
  
  # --- 你的靶点区 ---
  "Fpr1", "Fpr2",     # 细菌雷达
  "Ccr2", "Cxcr4"     # 迁移受体
)

p_specificity <- DotPlot(small_obj, features = check_list) + 
  RotatedAxis() +
  ggtitle("Target Specificity & Safety Check") +
  scale_color_gradient(low = "grey95", high = "red")

p_specificity
ggsave("01_project_GSE141259/04_output_plots/Specificity_Check.png", plot = p_specificity, width = 10, height = 6)

print("🎉 完美！这张图就是你要给 Tobias 看的‘安全性评估’！")



print("📊 正在进行统计学验证 (Finding Markers)...")

# 我们来算算 Cluster 1 (AM) 相比于所有其他细胞，哪些基因是显著高表达的？
# ident.1 = 1 (AM的编号，请根据你之前的图确认是1还是0)
am_markers <- FindAllMarkers(small_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 我们只看 Cluster 1 (AM) 的前 10 名大哥
top10_am <- am_markers %>%
  filter(cluster == "1") %>%  # 假设你的 AM 是 Cluster 1
  top_n(n = 10, wt = avg_log2FC)

print(top10_am)

# 把这个表格存下来，可以拿着数字说话！
write.csv(top10_am, "01_project_GSE141259/04_output_plots/AM_Top_Markers_Stats.csv")

print("✅ 统计学证据已到手！")





print("🔥 加班模式开启！正在全面验证 CCR 家族 (Ccr1-5) + FPR 家族...")

# 1. 定义全套名单 (Full List)
# -------------------------------------------------------------
# A. 身份 Marker (用来定位细胞类型)
identity_markers <- c("Emr1", "Marco",    # AM (你的主角)
                      "Col1a1",           # 成纤维 (干扰项)
                      "Pecam1")           # 血管 (干扰项)

# B. 目标受体 (CCR 全家桶 + FPR 双煞)
# 咱们按顺序排好，方便对比
targets <- c("Fpr1", "Fpr2",              # 细菌雷达
             "Ccr1", "Ccr2", "Ccr3", "Ccr4", "Ccr5") # 趋化受体全家桶

# 合并成一个列表
all_genes <- c(identity_markers, targets)

# 2. 绘制超级气泡图 (Super DotPlot)
# -------------------------------------------------------------
p_full_screen <- DotPlot(small_obj, features = all_genes) + 
  RotatedAxis() +
  ggtitle("Comprehensive Receptor Screen: FPRs & CCR1-5") +
  scale_color_gradient(low = "grey95", high = "red") # 红色越深越强

# 显示图片
p_full_screen

# 3. 保存这张终极战果
# -------------------------------------------------------------
ggsave("01_project_GSE141259/04_output_plots/Full_CCR_FPR_Screen.png", plot = p_full_screen, width = 11, height = 6)

print("🎉 搞定！CCR1-5 和 Fpr1/2 已全部列队完毕！快去看图！")


# 上述步骤整理后的版本
# ==============================================================================
# 0. 准备工作 (Setup)
# ==============================================================================
# 清空环境 (可选，保持干净)
rm(list = ls())
gc()

# 加载必要的包
library(Seurat)
library(tidyverse)

# ...然后下面再接我刚才发给你的那三段代码...
# ==============================================================================
# 1. 读取数据 (Read Data)
# ==============================================================================

# 使用 ReadMtx 读取 10x 格式的三件套文件
# mtx = 表达量矩阵(数字), cells = 细胞ID, features = 基因名
raw_counts <- ReadMtx(
  mtx = "01_project_GSE141259/00_data_raw/GSE141259_WholeLung_rawcounts.mtx.gz",
  cells = "01_project_GSE141259/00_data_raw/GSE141259_WholeLung_barcodes.txt.gz",
  features = "01_project_GSE141259/00_data_raw/GSE141259_WholeLung_genes.txt.gz",
  feature.column = 1  # 关键修正：告诉R基因名在第1列 (默认找第2列会报错)
)

# ==============================================================================
# 2. 创建对象与初筛 (Create Object & QC)
# ==============================================================================

# 把矩阵装进 Seurat 对象这个“收纳箱”里
strunz_obj <- CreateSeuratObject(
  counts = raw_counts,
  project = "Strunz_Lung",
  min.cells = 3,       # 过滤基因：如果一个基因在少于3个细胞里表达，视为噪音扔掉
  min.features = 200   # 过滤细胞：如果一个细胞测到的基因少于200个，视为垃圾/碎片扔掉
)

# 内存清理：把刚才那一大盆原始数据倒掉，只保留装箱后的对象
rm(raw_counts)
gc()

# 抽样 (为了练习时不卡顿，真实跑数据时可跳过这一步)
set.seed(123) # 设置随机种子，保证每次抽的人都一样
small_obj <- subset(strunz_obj, cells = sample(Cells(strunz_obj), 5000))
# ==============================================================================
# 3. 标准处理流程 (Standard Pipeline)
# ==============================================================================

small_obj <- small_obj %>%
  
  # A. 标准化 (Normalization)
  # 目的：消除测序深度的影响 (把大家拉到同一起跑线)
  # 原理：LogNormalize (总数归一化后取对数)
  NormalizeData() %>%
  
  # B. 找高变基因 (Feature Selection)
  # 目的：挑出最能区分细胞类型的2000个“特征基因” (如 Emr1, Col1a1)
  # 忽略那些所有细胞都一样的“管家基因”
  FindVariableFeatures() %>%
  
  # C. 归一化 (Scaling)
  # 目的：把数据转换成 Z-score (均值为0，方差为1)
  # 作用：防止高表达基因(如胶原蛋白)权重过大，掩盖了其他基因
  ScaleData() %>%
  
  # D. 降维 PCA (Linear Dimensionality Reduction)
  # 目的：把2000个基因的复杂关系，压缩成50个主成分(PC)
  # 作用：提炼核心特征，去除噪音
  RunPCA() %>%
  
  # E. 降维 UMAP (Non-linear Dimensionality Reduction)
  # 目的：把50个PC压缩成2D坐标(x,y)，为了画那张漂亮的散点图
  # dims = 1:15 意思是用前15个主成分来画图 (通常够用了)
  RunUMAP(dims = 1:15) %>%
  
  # F. 找邻居 (Neighbors)
  # 目的：在数学空间里算算谁和谁离得近，构建社交网络
  FindNeighbors(dims = 1:15) %>%
  
  # G. 聚类 (Clustering)
  # 目的：把抱团的邻居圈起来，贴上 0,1,2,3 的标签
  # resolution = 0.5 是分辨率：数值越大，分得越细(群越多)
  FindClusters(resolution = 0.5)
# ==============================================================================
# 4. 可视化 (Plotting) - 修正豪华版
# ==============================================================================

# A. 基础地图 (UMAP)
# 看看细胞分成了几大坨。留着它，用来对照看Cluster编号！
p1 <- DimPlot(small_obj, label = TRUE) + ggtitle("Strunz Lung (Cluster Map)")

# B. 终极气泡图 (Super DotPlot)
# 这次我们把名单列全！
features_list <- c(
  # --- 1. 身份区 (Identity) ---
  "Emr1", "Marco",    # AM (主角)
  "Col1a1",           # 成纤维 (干扰)
  "Pecam1",           # 血管 (干扰)
  
  # --- 2. 靶点区 (Targets) ---
  "Fpr1", "Fpr2",     # 细菌雷达
  "Ccr1", "Ccr2", "Ccr3", "Ccr4", "Ccr5" # CCR 全家桶 (上次漏掉的就在这！)
)

p_screen <- DotPlot(small_obj, features = features_list) + 
  RotatedAxis() + 
  scale_color_gradient(low = "grey95", high = "red") +
  ggtitle("Comprehensive Screen: Identity + FPRs + CCR1-5")

# 显示图片
print(p1)       # 第一张：地图
print(p_screen) # 第二张：你要的那个全家桶！

# 保存图片 (覆盖之前的旧图)
ggsave("01_project_GSE141259/04_output_plots/Final_Full_Screen.png", plot = p_screen, width = 12, height = 6)
ggsave("01_project_GSE141259/04_output_plots/Final_UMAP.png", plot = p1)


# 分步骤进行
print("🎨 2. 正在绘制单基因表达 UMAP (FeaturePlots) - 防弹版...")

# --- 2.1 身份鉴定组 (Identity) ---
identity_genes <- c("Emr1", "Marco", "Col1a1", "Pecam1")

# 修正：去掉了那个导致报错的 plot_annotation
p2_identity <- FeaturePlot(small_obj, features = identity_genes, ncol = 2) 

print(p2_identity)
ggsave("01_project_GSE141259/04_output_plots/2_Identity_FeaturePlots.png", plot = p2_identity, width = 10, height = 8)


# --- 2.2 靶点全家桶组 (Targets: Fpr + Ccr) ---
target_genes <- c("Fpr1", "Fpr2", 
                  "Ccr1", "Ccr2", "Ccr3", "Ccr4", "Ccr5")

# 修正：去掉了那个导致报错的 plot_annotation
p2_targets <- FeaturePlot(small_obj, features = target_genes, ncol = 3) 

print(p2_targets)
ggsave("01_project_GSE141259/04_output_plots/2_Targets_FeaturePlots.png", plot = p2_targets, width = 12, height = 12)

print("✅ 第二部分 FeaturePlot 搞定！绝不报错！")
print("📊 3. 正在绘制汇总点图 (DotPlot) - 最终统计篇...")

# ==============================================================================
# 定义基因列表 (全家桶)
# ==============================================================================
# 我们按逻辑排个序，这样画出来的图好看
all_genes_list <- c(
  # 1. 身份区 (定坐标)
  "Emr1", "Marco",    # AM 主角
  "Col1a1", "Pecam1", # 干扰项
  
  # 2. FPR 家族 (细菌雷达)
  "Fpr1", "Fpr2", 
  
  # 3. CCR 家族 (趋化受体全家桶)
  # 按数字顺序排好
  "Ccr1", "Ccr2", "Ccr3", "Ccr4", "Ccr5"
)

# ==============================================================================
# 绘图 (DotPlot)
# ==============================================================================
# scale_color_gradient: 设定颜色，grey95是浅灰(低表达)，red是红(高表达)
p3_dotplot <- DotPlot(small_obj, features = all_genes_list) + 
  RotatedAxis() + 
  scale_color_gradient(low = "grey95", high = "red") +
  ggtitle("Comprehensive DotPlot: Identity + FPRs + CCRs")

# 显示图片 (Review)
print(p3_dotplot)

# ==============================================================================
# 保存 (Save)
# ==============================================================================
ggsave("01_project_GSE141259/04_output_plots/3_Final_DotPlot.png", plot = p3_dotplot, width = 12, height = 6)

print("🎉 全部搞定！Map (地图), Blueprints (蓝图), Bubbles (气泡) 全部集齐！")
print("📂 快去 04_output_plots 文件夹检阅你的战利品！")



print("🕵️‍♀️ 开启盲盒模式！正在寻找 AM 的独家秘籍...")

# ==============================================================================
# 寻找 AM (Cluster 1) 的特异性 Marker
# ==============================================================================

# FindMarkers 是 Seurat 最强大的功能之一
# ident.1 = 1  ->  我们要找 Cluster 1 (AM)
# min.pct = 0.25 -> 只看那些至少在 25% 的 AM 里表达的基因 (过滤掉噪音)
# only.pos = TRUE -> 只看“高表达”的，不看“低表达”的

am_markers <- FindMarkers(small_obj, ident.1 = 1, min.pct = 0.25, only.pos = TRUE)

# ==============================================================================
# 整理排行榜 (Top 10)
# ==============================================================================
# 我们按 avg_log2FC (倍数变化) 排序，看看谁是那个“最靓的仔”
top10_am <- am_markers %>%
  arrange(desc(avg_log2FC)) %>% # 从高到低排
  head(10) # 只看前 10 名

print("🏆 AM 里的 Top 10 高表达基因是：")
print(top10_am)

# ==============================================================================
# 把前 10 名存下来，或者画个图看看
# ==============================================================================
# 我们把这些自动找出来的 Top 基因画个气泡图
# rownames(top10_am) 就是那 10 个基因的名字
p_discovery <- DotPlot(small_obj, features = rownames(top10_am)) + 
  RotatedAxis() +
  ggtitle("Top 10 Genes Defining AM (Unbiased Discovery)")

print(p_discovery)
ggsave("01_project_GSE141259/04_output_plots/AM_Top10_Discovery.png", plot = p_discovery, width = 12, height = 6)