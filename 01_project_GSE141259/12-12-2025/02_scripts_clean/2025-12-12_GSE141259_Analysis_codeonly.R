
rm(list = ls())
gc()

library(Seurat)
library(tidyverse)

# 1. 读取数据 (Read Data)
raw_counts <- ReadMtx(
  mtx = "01_project_GSE141259/00_data_raw/GSE141259_WholeLung_rawcounts.mtx.gz",
  cells = "01_project_GSE141259/00_data_raw/GSE141259_WholeLung_barcodes.txt.gz",
  features = "01_project_GSE141259/00_data_raw/GSE141259_WholeLung_genes.txt.gz",
  feature.column = 1  # 关键修正：告诉R基因名在第1列 (默认找第2列会报错)
)

# 📂 单细胞测序数据读取：三种常见模式总结
# 1）. 标准 10x Genomics 格式（文件夹模式）
data <- Read10X(data.dir = "路径/filtered_feature_bc_matrix/")
seu_obj <- CreateSeuratObject(counts = data)

# 2）. HDF5 格式（单文件模式）
data <- Read10X_h5(filename = "路径/data.h5")
seu_obj <- CreateSeuratObject(counts = data)

# 3）. GEO/自定义格式
data <- ReadMtx(
  mtx = "路径/GSExxxx_matrix.mtx.gz",      # 数值文件
  cells = "路径/GSExxxx_barcodes.txt.gz",  # 列名文件
  features = "路径/GSExxxx_genes.txt.gz",  # 行名文件
  feature.column = 1  # ⚠️关键点：如果 genes 文件只有一列，必须加这个参数！
)
# 注意：这是最灵活的方法，但也最容易因为参数设置错误（如列数不对）而报错。


# 2. 创建对象与初筛 (Create Object & QC)

strunz_obj <- CreateSeuratObject(
  counts = raw_counts,
  project = "Strunz_Lung",
  min.cells = 3,       # 过滤基因：如果一个基因在少于3个细胞里表达，视为噪音扔掉
  min.features = 200   # 过滤细胞：如果一个细胞测到的基因少于200个，视为垃圾/碎片扔掉
)

# 内存清理
rm(raw_counts)
gc()

# 抽样 (为了练习时不卡顿，真实跑数据时可跳过这一步)
set.seed(123) # 设置随机种子，保证每次抽的人都一样
small_obj <- subset(strunz_obj, cells = sample(Cells(strunz_obj), 5000))


# 3. 标准处理流程 (Standard Pipeline)

small_obj <- small_obj %>%
  
  NormalizeData() %>%

  FindVariableFeatures() %>%
  
  ScaleData() %>%
  
  RunPCA() %>%
  
  RunUMAP(dims = 1:15) %>%
  
  FindNeighbors(dims = 1:15) %>%
  
  FindClusters(resolution = 0.5)


# 4. 可视化 (Visualization)

# 图 A: 基础地图 (Cluster Map)
p1_clusters <- DimPlot(small_obj, label = TRUE, label.size = 5) + 
  ggtitle("1. Cluster Map (Who is where?)")

# 图 B: 身份鉴定图 (Identity FeaturePlots)
p2_identity <- FeaturePlot(small_obj, 
                           features = c("Emr1", "Marco",    # AM (主角)
                                        "Col1a1", "Pecam1"), # 干扰项
                            ncol = 2) # 排成 2列 x 2行

# 图 C: 靶点验证图 (Target FeaturePlots)
p3_targets <- FeaturePlot(small_obj, 
                          features = c("Fpr1", "Fpr2", "Ccr2", "Ccr5"),
                          ncol = 2,
                          # 增加一个小技巧: 调整颜色 (浅灰 -> 红)，红色可能比蓝色更显眼
                          cols = c("lightgrey", "red")) 

# 图 D: 那个超级气泡图 (Super DotPlot) - 统计学汇总
features_list <- c("Emr1", "Marco", "Col1a1", "Pecam1", # 身份
                   "Fpr1", "Fpr2", "Ccr1", "Ccr2", "Ccr5") # 靶点
p4_dotplot <- DotPlot(small_obj, features = features_list) + 
  RotatedAxis() + 
  scale_color_gradient(low = "grey95", high = "red") +
  ggtitle("4. Statistical Summary")

# 5. 保存战果 (Save)
# ==============================================================================
ggsave("01_project_GSE141259/04_output_plots/Map_Clusters.png", plot = p1_clusters)
ggsave("01_project_GSE141259/04_output_plots/Map_Identity_Spatial.png", plot = p2_identity, width = 10, height = 8)
ggsave("01_project_GSE141259/04_output_plots/Map_Targets_Spatial.png", plot = p3_targets, width = 10, height = 8)
ggsave("01_project_GSE141259/04_output_plots/Stats_DotPlot.png", plot = p4_dotplot, width = 12, height = 6)
# 练习时可以不用ggsave，直接在右下角export功能输出图片，避免了每次run all又重新保存一堆图片，若正式用于发文章，请用ggsave：
