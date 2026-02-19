# ==============================================================================
# Script: 13_NicheNet_Analysis_Final_Fix.R
# Purpose: 推断配体 (带 Mouse->Human 基因名转换)
# ==============================================================================

library(Seurat)
library(nichenetr)
library(tidyverse)
library(ggplot2)

# 1. 读取数据和知识库
# ------------------------------------------------------------------------------
rds_path <- "lianxi/01_data_processed/lung_obj_final_analysis.rds"
if(!exists("seurat_obj")) seurat_obj <- readRDS(rds_path)

message("📂 正在加载 NicheNet 知识库...")
ligand_target_matrix <- readRDS("lianxi/ligand_target_matrix.rds")
lr_network <- readRDS("lianxi/lr_network.rds")
weighted_networks <- readRDS("lianxi/weighted_networks.rds")

# 2. 设置正确的细胞身份
# ------------------------------------------------------------------------------
message("🔧 校正细胞身份...")
if("cell.type" %in% colnames(seurat_obj@meta.data)) {
  Idents(seurat_obj) <- "cell.type"
}
# 确认细胞类型
sender <- c("Fibroblasts", "Myofibroblasts")
sender <- sender[sender %in% unique(Idents(seurat_obj))] # 确保只选存在的
receiver <- "Krt8 ADI"

message("🎯 锁定目标：Sender [", paste(sender, collapse="+"), "] -> Receiver [", receiver, "]")

# 3. 提取基因并进行【物种转换】(Mouse -> Human)
# ------------------------------------------------------------------------------
message("🧬 正在提取基因并翻译 (Mouse -> Human)...")

# 定义一个简单的转换函数 (把首字母大写变成全大写)
mouse_to_human <- function(genes) {
  return(toupper(genes)) 
}

# 提取 Mouse 基因
genes_receiver_mouse <- get_expressed_genes(receiver, seurat_obj, pct = 0.10)
genes_sender_mouse <- get_expressed_genes(sender, seurat_obj, pct = 0.10)

# 翻译成 Human 基因 (NicheNet 只认识 Human)
expressed_genes_receiver <- mouse_to_human(genes_receiver_mouse)
expressed_genes_sender <- mouse_to_human(genes_sender_mouse)

# 确保翻译后的基因真的在数据库里
expressed_genes_receiver <- expressed_genes_receiver[expressed_genes_receiver %in% colnames(ligand_target_matrix)]
# 对于 sender，我们要看它是否在网络中作为 ligand 存在
all_ligands <- lr_network$from %>% unique()
expressed_genes_sender <- expressed_genes_sender[expressed_genes_sender %in% all_ligands]

# 4. 定义 "Genes of Interest" (ADI 的变化)
# ------------------------------------------------------------------------------
message("📝 提取 ADI 特征基因...")

# 找差异基因 (Mouse 数据)
deg_adi <- FindMarkers(seurat_obj, ident.1 = receiver, ident.2 = "AT2 cells", 
                       only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 翻译差异基因
genes_of_interest_mouse <- rownames(deg_adi)
genes_of_interest_human <- mouse_to_human(genes_of_interest_mouse)

# 过滤：必须是数据库里有的靶基因
geneset_human <- genes_of_interest_human[genes_of_interest_human %in% colnames(ligand_target_matrix)]

message("✅ 翻译前(Mouse): ", length(genes_of_interest_mouse), " -> 翻译后匹配(Human): ", length(geneset_human))

if(length(geneset_human) == 0) stop("❌ 还是没找到匹配基因！请检查网络或数据。")

# 5. 预测配体活动 (v2.0)
# ------------------------------------------------------------------------------
message("🚀 正在推理潜在配体...")

potential_ligands <- lr_network %>% 
  filter(from %in% expressed_genes_sender & to %in% expressed_genes_receiver) %>%
  pull(from) %>% unique()

# 这里的 geneset 参数用的是翻译后的 Human 基因
ligand_activities <- predict_ligand_activities(geneset = geneset_human,
                                               background_expressed_genes = expressed_genes_receiver,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

best_upstream_ligands <- ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

message("🏆 嫌疑人名单 (Top Ligands): ", paste(head(best_upstream_ligands, 5), collapse = ", "))

# 6. 可视化
# ------------------------------------------------------------------------------
message("🎨 正在绘图...")

# 图1: 配体排名
p_ligand_activity <- ligand_activities %>%
  top_n(20, pearson) %>%
  ggplot(aes(x = pearson, y = reorder(test_ligand, pearson))) + 
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(title = "Top Predicted Ligands (Sender: Fibroblasts)", x = "Pearson Correlation", y = "Ligand") +
  theme_classic()

ggsave("lianxi/04_output_plots/NicheNet_Ligand_Activity.png", p_ligand_activity, width = 6, height = 8)

# 图2: 调控热图
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links, geneset = geneset_human, ligand_target_matrix = ligand_target_matrix, n = 200) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25) #稍微放宽一点阈值

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links))
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()

p_ligand_target <- vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized Ligands", "ADI Target Genes", 
                      color = "purple", legend_position = "top", 
                      x_axis_position = "top", legend_title = "Regulatory Potential") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0)) + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01))

ggsave("lianxi/04_output_plots/NicheNet_Ligand_Target_Heatmap.png", p_ligand_target, width = 12, height = 8)

message("🎉 大功告成！快去看图！")


# ==============================================================================
# 📊 Figure Interpretation: NicheNet Ligand-Target Analysis (图解指南)
# ==============================================================================

# ------------------------------------------------------------------------------
# Figure 1: NicheNet_Ligand_Activity.png (Top Predicted Ligands)
# Focus: What is the Fibroblast niche secreting? (成纤维细胞在分泌什么？)
# ------------------------------------------------------------------------------

# 1. The TGF-beta Activator (纤维化主开关):
#    - [Observation]: THBS1 (Thrombospondin 1) is the top-ranked predicted ligand.
#    - [Meaning]: THBS1 is a major endogenous activator of latent TGF-beta. 
#      Fibroblasts are actively establishing a highly pro-fibrotic signaling environment.
#    - [中文]: THBS1 高居榜首。它是 TGF-β（纤维化核心因子）的强效激活剂，说明成纤维细胞正在主动制造纤维化风暴。

# 2. Extracellular Matrix Remodeling (基质重塑/僵硬化):
#    - [Observation]: Enrichment of structural ECM proteins (COL4A1, LAMB1, FBN1, FN1).
#    - [Meaning]: Validates the CellChat findings. The ADI cells are trapped in a stiff, 
#      pathological extracellular matrix, which likely acts as a mechanical stressor driving their reprogramming.

# ------------------------------------------------------------------------------
# Figure 2: NicheNet_Ligand_Target_Heatmap.png (Regulatory Potential)
# Focus: What do these ligands DO to the ADI cells? (这些配体导致了什么后果？)
# ------------------------------------------------------------------------------

# 1. The SPP1-ICAM1 Axis (炎症与衰老的纽带):
#    - [Observation]: Strong regulatory potential between Fibroblast-derived SPP1 (and LAMB1/THBS1) 
#      and the target gene ICAM1 in ADI cells.
#    - [Meaning]: ICAM1 is a critical adhesion molecule and a well-known marker of 
#      cellular senescence and the SASP (Senescence-Associated Secretory Phenotype).
#    - [Conclusion]: Signals from the fibrotic niche directly instruct ADI cells to adopt 
#      a pro-inflammatory, senescent phenotype (upregulating ICAM1), preventing their normal regeneration.
#    - [中文]: 这是一个关键发现。成纤维细胞的配体直接驱动了 ADI 细胞中 ICAM1 的表达。ICAM1 是经典的衰老/促炎标志物，这意味着微环境直接“锁死”了 ADI 的病理状态，迫使它们成为炎症推手。

# ==============================================================================
# 📝 Summary for Manuscript (论文结论):
# "NicheNet analysis reveals that the fibrotic niche (Fibroblasts) drives the ADI state 
#  through the secretion of ECM components (THBS1, LAMB1) and SPP1, which subsequently 
#  upregulate the senescence and adhesion marker ICAM1 on transitioning epithelial cells."
# ==============================================================================
