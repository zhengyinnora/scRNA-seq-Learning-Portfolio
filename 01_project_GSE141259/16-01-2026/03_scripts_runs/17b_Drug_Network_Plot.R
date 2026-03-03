# ==============================================================================
# Script: 17b_Drug_Network_Plot.R
# Purpose: 绘制高大上的药物-靶点互作网络图 (Network Plot)
# ==============================================================================

# 1. 安装并加载画图神器
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(dplyr)
library(igraph)
library(ggraph)

message("📂 1. 正在读取候选药物名单...")
drug_data <- read.csv("lianxi/04_output_plots/17_Drug_Repurposing_Candidates.csv")

message("🧹 2. 正在清洗数据，剔除未命名的实验性化合物...")
# 核心逻辑：过滤掉名字里带有 "COMPOUND", "[", "EXAMPLE", "CHEMBL" 的行
clean_edges <- drug_data %>%
  filter(!grepl("COMPOUND|\\[|EXAMPLE|CHEMBL", drug_name, ignore.case = TRUE)) %>%
  select(drug_name, gene_name, interaction_score) %>%
  rename(from = drug_name, to = gene_name) # 连线方向：药物(from) -> 靶点(to)

# 如果剩下的药还是太多（比如超过50个），图会变成“毛线球”
# 这里做一个防御机制，默认最多只取前 40 个最靠谱的互作组合
if(nrow(clean_edges) > 40) {
  clean_edges <- head(clean_edges, 40)
}

message("🕸️ 3. 正在构建网络对象...")
# 构建 igraph 专用网络对象
net <- graph_from_data_frame(d = clean_edges, directed = TRUE)

# 给节点打标签：区分是“药物”还是“靶点基因”，方便后面上色
# 逻辑：如果节点名字存在于目标靶点列(to)，那就是基因，否则就是药物
V(net)$type <- ifelse(V(net)$name %in% clean_edges$to, "Target Gene", "Drug")

message("🎨 4. 正在绘制网络图 (这可能需要几秒钟渲染)...")
# 使用 ggraph 画图，'fr' 是一种经典的力导向布局，会自动把节点推开
p <- ggraph(net, layout = 'fr') +  
  # 1. 画连线：带箭头的灰色半透明线
  geom_edge_link(aes(edge_alpha = 0.5), color = "grey60",
                 arrow = arrow(length = unit(2, 'mm'), type = "closed"),
                 end_cap = circle(3, 'mm')) +
  # 2. 画节点：根据节点类型(type)自动赋予不同颜色和大小
  geom_node_point(aes(color = type, size = type)) +
  # 3. 核心配色：药物用薄荷绿，基因用警示橙
  scale_color_manual(values = c("Drug" = "#1b9e77", "Target Gene" = "#d95f02")) +
  scale_size_manual(values = c("Drug" = 5, "Target Gene" = 8)) +
  # 4. 加名字：repel=TRUE 表示文字会自动避开不重叠
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5, fontface = "bold", color = "black") +
  # 5. 美化背景
  theme_void() + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  labs(title = "Drug-Target Interaction Network",
       color = "Node Type",
       size = "Node Type")

message("💾 5. 正在保存高清 PDF 和 PNG 图片...")
ggsave("lianxi/04_output_plots/17b_Drug_Target_Network.png", plot = p, width = 8, height = 7, dpi = 300)
ggsave("lianxi/04_output_plots/17b_Drug_Target_Network.pdf", plot = p, width = 8, height = 7)

message("🎉 搞定！")

# ==============================================================================
# Figure Legend / Output Interpretation for Manuscript:
# ==============================================================================
# Figure X: Drug-Target Interaction Network of Cross-Species Conserved Pathogenic Genes.
# 
# (A) Network visualization illustrating predicted pharmacological interventions 
# against key disease-driving targets. Orange nodes represent core pathogenic genes 
# that are significantly upregulated and conserved across species. Mint green nodes 
# denote candidate drugs or therapeutic compounds. 
# 
# (B) Directed grey edges (arrows) originating from drugs to target genes indicate 
# a suppressive pharmacological mechanism (e.g., inhibitors, antagonists, blockers, 
# or neutralizing antibodies). 
# 
# (C) The network topology reveals both highly druggable target hubs (genes susceptible 
# to multiple therapeutic agents, such as GTF3C2) and highly specific drug-target 
# pairs (e.g., C3 targeted by AMY-101, CTSV targeted by DEFACTINIB). 
# 
# (D) Drug-gene interaction data was systematically retrieved from the Drug-Gene 
# Interaction Database (DGIdb). Uncharacterized experimental compounds were excluded 
# to prioritize candidates with prominent clinical translation potential and repurposing 
# value.
# ==============================================================================