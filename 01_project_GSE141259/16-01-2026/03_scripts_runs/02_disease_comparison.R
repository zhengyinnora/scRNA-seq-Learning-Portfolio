# ==============================================================================
# Script: 02_disease_comparison.R
# [CN] 这是你的第二个脚本，专门用来做分析。
#      第一步必须是“复活”上一关存好的数据。
# ==============================================================================

# 1. 还是要先加载这几个老伙计 (因为新脚本里环境是空的)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# 2. 告诉 R，你的存档文件存在哪个文件夹里
# (记得吗？刚才咱们定义过这个路径，新脚本里要再定义一次)
input_dir <- "lianxi/01_data_processed"

# 3. 关键一步：读取 .rds 文件！
# [EN] Load the saved object. 
# [CN] 就像游戏读档一样。
#      注意：括号里的文件名必须和你刚才保存的文件名一模一样！
#      readRDS 读取进来只是数据，我们需要用 "<-" 把这个数据赋值给 lung_obj 变量
lung_obj <- readRDS(file.path(input_dir, "lung_obj_processed_with_celltype.rds"))

# 4. 检查一下是不是复活成功了
# 如果输出显示 "An object of class Seurat..." 还有 29297 个细胞，那就是成功了！
print(lung_obj)

# ==============================================================================
# Step 1: Check Experimental Groups
# [CN] 检查分组情况。
#      我们使用 table() 函数来统计 'grouping' 这一列里到底有哪些组，以及每组有多少细胞。
# [EN] Inspect the experimental groups.
#      Use table() to count the number of cells in each group within the 'grouping' column.
# ==============================================================================

# 看看 grouping 这一列里都有啥
table(lung_obj$grouping)

# ==============================================================================
# Step 2: Plot the Time Course
# [CN] 既然发现了7个时间点，我们就把它们全画出来。
#      split.by = "grouping": 按照 d3, d7... PBS 自动拆分。
#      ncol = 3: 每行放3张图，排版好看点。
# ==============================================================================

# 1. 设置保存路径 (还是那个 plots 文件夹)
plots_dir <- "lianxi/04_output_plots"

# 2. 画图！(这步会生成7张小图拼在一起)
p_split <- DimPlot(lung_obj, reduction = "umap", group.by = "cell.type", 
                   split.by = "grouping", label = TRUE, repel = TRUE, ncol = 3) + NoLegend()

# 3. 保存高清大图 (因为有7张图，我们要把画布设得很大)
ggsave(file.path(plots_dir, "lianxixi_TimeCourse_Split.png"), plot = p_split, width = 20, height = 15)

# 4. 在右下角预览一下
p_split

# ==============================================================================
# Step 3: Reorder Time Points
# [CN] 修正时间顺序。
#      R 默认按字母排序 (d1, d10, d2...)，我们需要强制指定符合生物学逻辑的顺序。
#      Levels: PBS -> d3 -> d7 -> d10 -> d14 -> d21 -> d28
# ==============================================================================

# 1. 定义正确的顺序列表
time_order <- c("PBS", "d3", "d7", "d10", "d14", "d21", "d28")

# 2. 使用 factor 函数重塑 grouping 列
# 这步操作不会改变数据，只会改变它们的“排队顺序”
lung_obj$grouping <- factor(lung_obj$grouping, levels = time_order)

# ==============================================================================
# Visualization (Preview Mode)
# [CN] 再次画图预览。
#      这次你会发现小图的排列顺序变顺了！
#      注意：这次不自动保存，只在右下角显示。
# ==============================================================================

DimPlot(lung_obj, reduction = "umap", group.by = "cell.type", 
        split.by = "grouping", label = TRUE, repel = TRUE, ncol = 4) + NoLegend()

# [Optional] 如果觉得好看，你再手动运行下面这行来保存：
ggsave(file.path(plots_dir, "lianxixi_TimeCourse_Ordered.png"), width = 20, height = 10)

# ==============================================================================
# Step 4: Cell Type Proportion Analysis
# [CN] 细胞比例分析。
#      UMAP 只能看“位置”，柱状图才能看“数量”。
#      我们要计算每种细胞在不同天数占的百分比。
# ==============================================================================

# 1. 只是画图预览 (Preview)
# fill = "cell.type": 用颜色区分细胞类型
# x = "grouping": 横坐标是时间点
# position = "fill": 堆叠显示（总高100%）
ggplot(lung_obj@meta.data, aes(x = grouping, fill = cell.type)) +
  geom_bar(position = "fill") +
  theme_classic() +
  labs(y = "Proportion", x = "Time Point", title = "Cell Type Changes over Time") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # 让横坐标文字斜着放，防重叠

# [Optional] 觉得好看再存：
ggsave(file.path(plots_dir, "lianxixi_BarPlot_Proportion.png"), width = 12, height = 8)

# ==============================================================================
# Step 5: Cell Type Trend Lines (The "Stock Market" View)
# [CN] 细胞比例折线图。
#      既然柱状图颜色太乱，我们把感兴趣的“明星细胞”拎出来，画出它们随时间变化的趋势。
#      这就像看股票走势一样清晰。
# [EN] Generate line plots for specific cell type proportions over time.
# ==============================================================================

# 1. 先计算出每个时间点、每种细胞的比例表
prop_table <- lung_obj@meta.data %>%
  group_by(grouping, cell.type) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

# 2. 挑选几个我们最想看的“主角” (根据刚才的读图)
# 你可以把名字换成任何你想看的细胞
targets <- c("AM (PBS)", "AM (Bleo)", "Krt8 ADI", "Fibroblasts", "AT2 cells")

# 3. 只筛选这些主角来画图
prop_filtered <- prop_table %>% filter(cell.type %in% targets)

# 4. 画折线图！
p_line <- ggplot(prop_filtered, aes(x = grouping, y = proportion, color = cell.type, group = cell.type)) +
  geom_line(size = 1.5) +  # 线画粗点
  geom_point(size = 3) +   # 点画大点
  theme_classic() +
  labs(title = "Major Cell Type Dynamics", y = "Proportion of Total Cells")

# 5. 预览
p_line

# [Optional] 保存
ggsave(file.path(plots_dir, "lianxixi_LinePlot_Trends.png"), plot = p_line, width = 8, height = 6)

# ==============================================================================
# Step 6: Feature Plot (Gene Expression Heatmap)
# [CN] 基因特征图。
#      刚才我们看的是“细胞在哪”，现在我们要看“基因在哪”。
#      我们要验证：
#      1. Sftpc (AT2 标志物): 应该在健康的 AT2 细胞里高表达。
#      2. Krt8 (损伤标志物): 应该在那群新出现的 Krt8 ADI 细胞里高表达。
# [EN] Visualize expression of specific genes on the UMAP.
# ==============================================================================

# 1. 定义我们要查的基因 (这两个是经典的肺部 Marker)
my_features <- c("Sftpc", "Krt8")

# ------------------------------------------------------------------------------
# Part A: 先看“精华版” (只看 PBS vs d14) —— 预览方便！
# ------------------------------------------------------------------------------
# 创建一个小对象用于快速预览
small_obj <- subset(lung_obj, subset = grouping %in% c("PBS", "d14"))

# [Action 1: 预览] 画图！(红色的点代表基因高表达)
p_small <- FeaturePlot(small_obj, features = my_features, split.by = "grouping", 
                       order = TRUE, cols = c("lightgrey", "red"))
print(p_small)  # <--- 这次先让你看！

# [Action 2: 保存] 看满意了吗？满意了运行这一行存下来：
ggsave(file.path(plots_dir, "lianxixi_FeaturePlot_Genes_Small.png"), plot = p_small, width = 12, height = 6)


# ------------------------------------------------------------------------------
# Part B: 再看“完整版” (所有 7 个时间点) —— 信息量大！
# ------------------------------------------------------------------------------
# [Action 3: 预览] 生成包含所有组的大图
# 注意：这图会很大，右下角可能会有点挤，但你能看到全貌
p_big <- FeaturePlot(lung_obj, features = my_features, split.by = "grouping", 
                     order = TRUE, cols = c("lightgrey", "red"), ncol = 4)
print(p_big)   # <--- 先看图！

# [Action 4: 保存] 同样，看满意了再存：
ggsave(file.path(plots_dir, "lianxixi_FeaturePlot_Genes_AllTime.png"), plot = p_big, width = 20, height = 10)


# ==============================================================================
# Step 7 (Fixed): Refined Violin Plot
# [CN] 优化小提琴图：只看“重点嫌疑人”。
#      之前细胞太多挤成一条线了。我们现在只筛选出 5 种上皮细胞来对比。
# ==============================================================================

# 1. 定义我们要看的“VIP 名单” (必须和图里的名字一模一样)
vip_cells <- c("AT2 cells", "Krt8 ADI", "AT1 cells", "Club cells", "Ciliated cells")

# 2. 从大对象里把这些细胞“提取”出来，生成一个小对象
# subset = cell.type %in% vip_cells: 意思是只保留 cell.type 在名单里的细胞
epithelial_obj <- subset(lung_obj, subset = cell.type %in% vip_cells)

# 3. 画图预览！(Preview)
# 只有5个柱子，这次肯定宽敞！
p_vln_clean <- VlnPlot(epithelial_obj, features = "Krt8", group.by = "cell.type", pt.size = 0) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + # 名字横着放，不用歪脖子看了
  NoLegend() +
  ggtitle("Krt8 Expression in Epithelial Cells")

# 4. 在右下角看图
print(p_vln_clean)

# ==============================================================================
# [Action] 觉得好看再保存：
# ==============================================================================
ggsave(file.path(plots_dir, "lianxixi_ViolinPlot_Krt8_Clean.png"), plot = p_vln_clean, width = 8, height = 6)

# ==============================================================================
# Final Save: Analysis Complete
# [CN] 保存最终分析结果。
#      这个 .rds 文件现在包含了我们所有的发现（分组、细胞类型、基因分析）。
# ==============================================================================

saveRDS(lung_obj, file = file.path(output_dir, "lung_obj_final_analysis.rds"))

message("恭喜 Nora！所有的分析成果都已安全保存！🎉")
