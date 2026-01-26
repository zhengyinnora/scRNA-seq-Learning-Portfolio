# ==============================================================================
# Script: 05_Trajectory_Analysis_Complete.R
# Purpose: Seurat -> Monocle3 -> Trajectory -> Pseudotime (Full Workflow)
# ==============================================================================

# 1. 加载必要的包
library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(igraph) # 用来处理节点名字

# ==============================================================================
# Step 1: 数据准备 (Data Prep)
# ==============================================================================

# 读取 Seurat 对象
lung_obj <- readRDS("lianxi/01_data_processed/lung_obj_final_analysis.rds")
message("✅ 数据加载完成")

# 提取数据构建 CDS
data <- GetAssayData(lung_obj, assay = "RNA", slot = "counts")
cell_metadata <- lung_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(data))
rownames(gene_metadata) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

# ==============================================================================
# Step 2: 预处理 (Preprocessing & Transfer)
# ==============================================================================

# 1. 搬运 UMAP 坐标
reducedDims(cds)$UMAP <- Embeddings(lung_obj, reduction = "umap")

# 2. 【关键】搬运正确的细胞类型名字 (修复之前的 bug)
#    直接用 Seurat 里名为 "cell.type" 的那一列
colData(cds)$cell_type <- lung_obj@meta.data$cell.type

# 3. 运行聚类 (为了不报错)
message("⏳ 正在初始化聚类...")
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 4. 学习轨迹 (画出黑线)
message("⏳ 正在学习轨迹 (Learn Graph)...")
cds <- learn_graph(cds, use_partition = FALSE)

# ==============================================================================
# Step 3: 寻找起点 & 计算拟时序 (The Logic that Worked)
# ==============================================================================

message("📍 正在寻找起点 (Root Node)...")

# A. 找到所有 "AT2 cells" 所在的节点位置
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
at2_cell_ids <- colnames(cds)[colData(cds)$cell_type == "AT2 cells"]
at2_nodes <- closest_vertex[at2_cell_ids, ]

# B. 找到细胞最多的那个节点 (得到的是数字 ID，比如 "274")
root_node_id <- names(which.max(table(at2_nodes)))

# C. 【关键修复】把数字 ID 转换成 Monocle 真正的节点名 (比如 "Y_274")
graph_node_names <- igraph::V(cds@principal_graph[["UMAP"]])$name
root_node_name <- graph_node_names[as.numeric(root_node_id)]

message("✅ 找到起点！ID: ", root_node_id, " -> Name: ", root_node_name)

# D. 计算拟时序 (Pseudotime)
cds <- order_cells(cds, root_pr_nodes = root_node_name)

# ==============================================================================
# Step 4: 保存与画图 (Save & Plot)
# ==============================================================================

# 1. 保存这张最珍贵的图
p_final <- plot_cells(cds,
                      color_cells_by = "pseudotime",
                      label_cell_groups = FALSE,
                      label_leaves = FALSE,
                      label_branch_points = FALSE,
                      graph_label_size = 3) +
  ggtitle("Pseudotime Trajectory (Purple: Start -> Yellow: End)")

ggsave("lianxi/04_output_plots/Monocle3_Pseudotime_Final.png", p_final, width = 8, height = 6)
print(p_final)

# 2. 保存最终的 CDS 对象 (里程碑！)
#    以后要用直接读这个文件，不用再跑上面的代码了
saveRDS(cds, "lianxi/01_data_processed/lung_monocle_final.rds")

message("🏆 全部大功告成！Nora 简直是生信战神！")

# ==============================================================================
# [Biological Interpretation / Figure Legend Draft]
# ==============================================================================
# 1. Trajectory Construction: 
#    Monocle3 successfully inferred a continuous trajectory spanning from 
#    naive AT2 cells to the injury-associated Krt8+ ADI state.
#
# 2. Pseudotime Analysis:
#    - Root (Pseudotime = 0): Defined as 'AT2 cells' (Dark Purple region).
#    - Terminus (High Pseudotime): Corresponds to 'Krt8 ADI' cells (Yellow region).
#
# 3. Key Finding:
#    The pseudotime gradient visualizes the progressive transdifferentiation 
#    of AT2 cells. The continuous path suggests a gradual loss of AT2 identity 
#    and acquisition of the ADI phenotype, rather than a discrete jump.
# ==============================================================================
# ==============================================================================
# [分析笔记]
# 1. 轨迹方向：成功构建了从正常 AT2 到 Krt8 ADI 的发育轨迹。
# 2. 拟时序含义：
#    - 紫色区域 (Pseudotime low) = 起始状态 (AT2)，代表未受损/稳态。
#    - 黄色区域 (Pseudotime high) = 终末状态 (Krt8 ADI)，代表损伤后的转分化结果。
# 3. 结论：
#    图示清晰地展示了 AT2 细胞存在向 ADI 细胞转化的“可塑性” (Plasticity)。
#    细胞不仅是分成了两堆，而是展示出了中间过渡的过程。
# ==============================================================================

# ==============================================================================
# Part 3: 寻找随拟时序变化的基因 (Finding Trajectory-dependent Genes)
# ==============================================================================

message("🕵️‍♀️ 正在侦测关键基因... (这一步可能有点慢，请喝口水 ☕️)")

# 1. 核心计算：graph_test (Moran's I 检验)
#    这个函数会找出那些在空间/轨迹上表达有规律的基因
#    neighbor_graph = "principal_graph" 表示我们要沿着画的那条黑线找
pr_graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# 2. 筛选显著的基因 (q_value < 0.05)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

message("✅ 计算完成！找到了 ", length(pr_deg_ids), " 个相关基因！")

# ==============================================================================
# Part 4: 看看 Top 基因是谁？ (Visualization)
# ==============================================================================

# 1. 按相关性 (morans_I) 排序，找出变化最剧烈的 Top 4 基因
#    morans_I 越高，说明这个基因沿着轨迹的变化越明显
top_genes <- pr_graph_test_res %>% 
  arrange(desc(morans_I)) %>% 
  head(4) %>% 
  rownames()

message("🏆 变化最剧烈的 Top 4 基因是：", paste(top_genes, collapse = ", "))

# 2. 画出这 4 个基因在轨迹上的表达图
p3 <- plot_cells(cds, 
                 genes = top_genes,
                 show_trajectory_graph = FALSE,
                 label_cell_groups = FALSE,
                 label_leaves = FALSE) +
  scale_color_viridis_c(option = "magma") + # 用岩浆色配色，看着比较高级
  ggtitle("Top Trajectory-Dependent Genes")

print(p3)
ggsave("lianxi/04_output_plots/Monocle3_Top_Genes.png", p3, width = 10, height = 8)

# 3. 再画一个具体基因看看 (比如 Krt8)
#    看看是不是真的在终点高表达？
p4 <- plot_cells(cds, 
                 genes = "Krt8",
                 show_trajectory_graph = TRUE,
                 label_cell_groups = FALSE,
                 label_leaves = FALSE) +
  ggtitle("Krt8 Expression along Trajectory")

print(p4)
ggsave("lianxi/04_output_plots/Monocle3_Gene_Krt8.png", p4, width = 8, height = 6)

message("🎉 基因分析完成！快看看图里是不是有的基因只在终点亮？")

# ==============================================================================
# Part 4: 基因模块分析 & 热图 (完整修复版)
# ==============================================================================

message("🕵️‍♀️ 正在筛选显著变化的基因...")

# 1. 过滤：只保留那些“显著变化”的基因 (q_value < 0.05)
#    (pr_graph_test_res 是刚才进度条跑出来的结果，必须存在内存里)
deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

message("✅ 找到显著变化的基因数: ", length(deg_ids))

# 2. 【关键补丁】计算 PCA
#    (Monocle3 分组需要这一步，之前报错就是因为缺这个)
message("🛠️ 正在打补丁：计算 PCA...")
cds <- preprocess_cds(cds, num_dim = 50)

# 3. 基因聚类 (Find Gene Modules)
#    把几千个基因分成几个“战队”
message("🧩 正在把基因分组 (Finding Modules)...")
gene_module_df <- find_gene_modules(cds[deg_ids,], resolution = 1e-2)

# 4. 绘制热图 (Heatmap)
message("🎨 正在绘制热图...")

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cell_type)

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

p_heatmap <- pheatmap::pheatmap(agg_mat,
                                scale="column", 
                                clustering_method="ward.D2",
                                main = "Gene Modules along Trajectory")

print(p_heatmap)
ggsave("lianxi/04_output_plots/Monocle3_Heatmap_Modules.png", p_heatmap, width = 10, height = 10)

# ==============================================================================
# Part 5: 抓取关键基因 (Who is with Krt8?)
# ==============================================================================

message("🕵️‍♀️ 正在侦查 Krt8 的同伙...")

# 1. 自动找到 Krt8 所在的那个 Module
krt8_module <- gene_module_df %>% 
  filter(id == "Krt8") %>% 
  pull(module)

message("👉 发现 Krt8 藏在 Module ", krt8_module, " 里！")

# 2. 把这个 Module 里的所有基因找出来
krt8_partners <- gene_module_df %>% 
  filter(module == krt8_module) %>% 
  pull(id)

# 3. 打印前 20 个看看
message("⚠️ 这个战队里共有 ", length(krt8_partners), " 个基因。前 20 个是：")
print(head(krt8_partners, 20))

# 4. 保存成表格 (这就是你发文章要用的基因列表！)
write.csv(krt8_partners, "lianxi/04_output_plots/Krt8_ADI_Module_Genes.csv")

message("🏆 全部搞定！Nora 可以下班了！🥂")

# ==============================================================================
# [Part 4 & 5 Interpretation: Gene Module Analysis]
# ==============================================================================
# 1. Analysis Logic (分析逻辑):
#    - Instead of analyzing genes individually, we grouped them into "Co-expression Modules"
#      using the 'find_gene_modules' algorithm.
#    - Genes within the same module share similar expression patterns along the trajectory,
#      suggesting they are co-regulated or functionally related.
#
# 2. Heatmap Visualization (热图解读):
#    - The heatmap (Part 4) visualizes how these gene modules turn on/off over pseudotime.
#    - Key Observation: Specific modules show distinct activation patterns at the 
#      Terminus (Krt8 ADI state), representing the core gene signature of the transdifferentiation.
#
# 3. Target Identification (核心发现):
#    - We identified the specific module containing the marker gene 'Krt8'.
#    - The genes in this list (saved in CSV) are co-expressed with Krt8.
#    - Biologically, these genes likely represent the molecular machinery driving 
#      the AT2-to-ADI transition (e.g., Krt19, Lgals3).
# ==============================================================================

# ==============================================================================
# [Nora 的人话笔记] 关于 Gene Module (基因模块)
# ==============================================================================
# 1. 这是在干嘛？
#    把 20,000 个基因通过“行为模式”进行分类。
#    就好比把全校学生按“社团”分组：足球社、合唱团、文学社...
#
# 2. 为什么要分 Module？
#    因为我们要找“坏人团伙”！
#    我们知道 Krt8 是坏人 (ADI marker)，但不知道它的同伙是谁。
#    通过分 Module，我们找到了和 Krt8 在同一个社团的所有基因。
#
# 3. 结果怎么看？
#    - CSV 表格里的名单 = Krt8 所在的那个“社团”的全员名单。
#    - 表格左边的数字 = 学号/序号 (不是社团名！)。
#    - 这份名单里的基因，就是我们要找的“核心突变基因群”！
# ==============================================================================