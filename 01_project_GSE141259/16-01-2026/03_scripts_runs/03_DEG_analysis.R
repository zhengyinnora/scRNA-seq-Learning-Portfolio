# ==============================================================================
# Script: 03_DEG_analysis.R
# Purpose: Find Differential Expressed Genes (DEGs) to understand molecular mechanisms.
# ==============================================================================
# ==============================================================================
# Step 1:
# 1. 加载包
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 2. 设置路径 (关键修改在这里！)
# 告诉 R：先进 lianxi 文件夹，再找 01_data_processed
input_dir  <- "lianxi/01_data_processed"
# 图片输出也顺便改对
plots_dir  <- "lianxi/04_output_plots" 

# 3. 读取存档
# 拼接路径 -> 读取文件
target_file <- file.path(input_dir, "lung_obj_final_analysis.rds")
message("正在读取文件: ", target_file)

lung_obj <- readRDS(target_file)

# 4. 检查一下
message("🎉 数据加载成功！")
print(lung_obj)

# ==============================================================================
# Step 2: Set Identity (Lock on Targets)
# [CN] 锁定目标。
#      我们要对比细胞类型，所以必须把身份(Identity)切换到 cell.type。
# ==============================================================================

# 1. 切换身份
Idents(lung_obj) <- "cell.type"

# 2. 检查名单 (Roll Call)
# 看看 "Krt8 ADI" 和 "AT2 cells" 是不是在名单里
# 这一步是为了防止名字写错导致后面报错
print("当前所有的细胞类型：")
print(levels(lung_obj))

# 3. 既然我们要存结果，先把存放表格的文件夹建好
# (在 lianxi 文件夹下新建一个 05_results_tables)
table_dir <- "lianxi/05_results_tables"
if(!dir.exists(table_dir)) dir.create(table_dir)

# 注释：简单说就是：
# 如果你没有这个文件夹 -> 它可以帮你自动新建一个，防止后面保存报错。
# 如果你已经有了（就像你现在的情况）-> 它就会直接跳过这一步，什么都不做，绝对不会覆盖或删除你里面已有的东西。
# 1. 定义地址：我们要把结果存在哪里？
# table_dir <- "lianxi/05_results_tables"
# 2. 智能判断：
# dir.exists(table_dir) -> 问电脑："这个文件夹存在吗？"
# ! -> 意思是 "不" (Not)
# 也就是："如果 (文件夹不！存在) ..."
# if(!dir.exists(table_dir)) 
# 3. 执行动作：
# "... 那就创建一个！"
# dir.create(table_dir)

# ==============================================================================
# Step 3: Find Differential Expressed Genes (The Battle)
# [CN] 寻找差异基因。
#      我们要对比 "Krt8 ADI" (生病组) 和 "AT2 cells" (健康组)。
#      ident.1 = "Krt8 ADI": 我们关注的主角（它高表达什么？）
#      ident.2 = "AT2 cells": 对照组（谁变过来的？）
# ==============================================================================

message("开始计算差异基因... (可能需要跑个几十秒，喝口水 ☕️)")

# 1. 运行 FindMarkers
# min.pct = 0.25: 只看那些至少在 25% 的细胞里表达的基因（剔除杂音）
# logfc.threshold = 0.25: 只看那些表达量变化比较明显的基因
deg_table <- FindMarkers(lung_obj, ident.1 = "Krt8 ADI", ident.2 = "AT2 cells", 
                         min.pct = 0.25, logfc.threshold = 0.25)

# 2. 看看战果 (前几名是谁？)
print("差异基因前 10 名：")
print(head(deg_table, 10))

# 3. 保存这张珍贵的清单
# write.csv: 把它存成 Excel 能打开的表格
write.csv(deg_table, file = file.path(table_dir, "DEG_Krt8_vs_AT2.csv"))

message("🎉 计算完成！清单已保存到 05_results_tables 文件夹！")
# ==============================================================================
# Step 3.1: Where is Krt8? (The Search)
# [CN] 专门查查 Krt8 排在哪里，数值是多少。
# ==============================================================================

# 1. 直接从大表格里提取 Krt8 这一行
krt8_stat <- deg_table["Krt8", ]
print("Krt8 的考试成绩：")
print(krt8_stat)

# 2. 顺便看看排在它前面的“正数大哥”都有谁（只看 Log2FC > 1 的）
# 这一步能帮你发现除了 Krt8 以外，还有谁是这群细胞的标志物
top_up_genes <- deg_table %>% 
  filter(avg_log2FC > 1) %>%  # 只看升高的
  arrange(desc(avg_log2FC))   # 从高到低排

print("升高基因里的前 10 名：")
print(head(top_up_genes, 10))

# ==============================================================================
# Step 4: Volcano Plot (Visual Confirmation)
# [CN] 这一步只画图，把 Krt8 标出来，让你亲眼看到它的位置。
# ==============================================================================

# 1. 挑选几个我们要点名的“明星基因”
# 把 Krt8 放第一个，必须要看到它！
# ==============================================================================
# [注释] 为什么选这 5 个基因做标记？(Selection Logic)
# ------------------------------------------------------------------------------
# 1. 身份互换证据 (Identity Switch):
#    - Krt8 (应上调): 这是这群细胞的“新身份证”，必须高表达。
#    - Sftpc (应下调): 这是健康 AT2 的“旧制服”，丢失说明去分化 (Dedifferentiation)。
#
# 2. 状态机制证据 (Mechanism):
#    - Cdkn1a (p21): 细胞周期“刹车片”。证明这群细胞处于损伤后的“停滞/衰老状态” (Senescence)。
#
# 3. 榜单高分证据 (Top Hits):
#    - S100a6 & Clu: 刚才榜单里的前十名“大哥”。代表细胞正在经历剧烈的压力应激和存活反应。
# ==============================================================================

# 1. 挑选几个我们要点名的“明星基因”
label_genes <- c("Krt8", "Sftpc", "Cdkn1a", "S100a6", "Clu")

# 2. 给表格加个标签
# 如果基因名字在名单里，就标出来，否则留空
deg_table$label <- ifelse(rownames(deg_table) %in% label_genes, rownames(deg_table), "")

# 3. 画火山图
p_volcano <- ggplot(deg_table, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = avg_log2FC > 0), alpha = 0.5, size = 1) + 
  geom_text(aes(label = label), vjust = -1, color = "black", fontface = "bold") + 
  scale_color_manual(values = c("blue", "red"), labels = c("AT2特征 (下调)", "Krt8 ADI特征 (上调)")) +
  theme_minimal() +
  labs(title = "Krt8 ADI vs AT2: 基因大乱斗")

# 4. 预览
print(p_volcano)

# 5. 保存
ggsave(file.path(plots_dir, "lianxixi_Volcano_Krt8_Final.png"), p_volcano, width = 8, height = 6)
# ==============================================================================
# [结果解读] 火山图真相：一场惨烈的身份互换与挣扎
# (Volcano Plot Interpretation: Identity Switch & Stress Struggle)
# ------------------------------------------------------------------------------
# 这张图清晰地展示了从健康 AT2 细胞到病态 Krt8 ADI 细胞的剧变过程：
#
# 1. 【旧身份的崩塌】(左上角蓝点):
#    - Sftpc: 作为 AT2 细胞最核心的标志物，它极显著地在 Krt8 ADI 中消失（深蓝色、高位置）。
#    -> 结论：细胞已经“去分化”，彻底放弃了制造表面活性物质的职能。
#
# 2. 【新身份的确立】(右侧中红点):
#    - Krt8 & Clu: 显著上调，位于右侧红色区域。
#    -> 结论：成功抓到了嫌疑人，确认这群细胞换上了新的损伤标志“马甲”。
#
# 3. 【痛苦的现状】(右上角顶端红点):
#    - S100a6 & Cdkn1a (p21): 这两个基因飞得最高，代表变化最剧烈、最显著。
#    -> 结论：这揭示了细胞的真实处境——它们正处于极度的“压力应激”(S100a6) 和“周期停滞”(Cdkn1a) 状态。它们虽然活着，但被锁死在了受损状态，无法修复。
# ==============================================================================

# ==============================================================================
# Step 5: Feature Plot (Visual Inspection on Map)
# [CN] 空间验证。
#      把刚才那几个明星基因，投射到 UMAP 地图上。
#      看看 Krt8 是不是真的只亮在那一小撮细胞上？
# ==============================================================================

# 1. 设置我们要看的基因列表
# Sftpc (AT2标志), Krt8 (损伤标志), Cdkn1a (停滞标志), S100a6 (压力标志)
features_to_plot <- c("Sftpc", "Krt8", "Cdkn1a", "S100a6")

# 2. 画图 (FeaturePlot)
# cols = c("lightgrey", "red"): 没表达是灰色，高表达是红色
p_feature <- FeaturePlot(lung_obj, features = features_to_plot, 
                         cols = c("lightgrey", "red"), ncol = 2)

# 3. 预览
print(p_feature)

# 4. 保存
ggsave(file.path(plots_dir, "lianxixi_FeaturePlot_KeyGenes.png"), p_feature, width = 10, height = 8)
# ==============================================================================
# [结果解读] UMAP 空间验证：嫌疑人的“藏身之处”
# (Step 5 Interpretation: Spatial Colocalization)
# ------------------------------------------------------------------------------
# 这张四宫格图 (FeaturePlot) 提供了无可辩驳的空间证据：
#
# 1. 【位置互斥】(Exclusivity):
#    - Sftpc (左上) 在大部分区域高表达，唯独在上方中间的一小撮细胞群中“熄灭”了。
#    - Krt8 (右上) 恰恰就在 Sftpc 熄灭的同一个位置“点亮”了。
#    -> 结论：这完美展示了 AT2 -> Krt8 ADI 的身份转换过程。
#
# 2. 【状态共定位】(Colocalization):
#    - Cdkn1a (左下) 和 S100a6 (右下) 的高亮区域，与 Krt8 的分布高度重合。
#    -> 结论：这证明了 Krt8 ADI 细胞群就是那群处于“衰老停滞”和“高压应激”状态的病理细胞。
#    -> 生物学意义：我们定位到了肺纤维化中一群具体的、功能失调的细胞亚群。
# ==============================================================================