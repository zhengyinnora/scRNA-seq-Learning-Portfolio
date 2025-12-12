# ==========================================
# Nora的单细胞练习脚本 (Version 2.0)
# 时间：2025-12-11
# ==========================================

# 1. 加载必要的包 (准备武器)
library(Seurat)
library(tidyverse)

# 2. 加载数据 (现在是练习数据，以后换成 nora_obj <- readRDS("师兄的数据.rds"))
nora_obj <- pbmc_small

# ----------------- 数据检查 (QC) -----------------

# 3. 看一下对象结构 (有多少细胞？多少基因？)
nora_obj

# 4. 看有哪些基因 (相当于看菜单，找找有没有感兴趣的)
# 注意：真实数据基因有几万个，屏幕会刷屏，可以用 head(rownames(nora_obj)) 只看前几个
rownames(nora_obj)

# ----------------- 可视化 (Plotting) -----------------

# 5. 按照 Cluster 画全景图 (UMAP/tSNE)
DimPlot(nora_obj)

# 6. 按照 Cluster 画 FeaturePlot (看特定基因在哪里表达)
# 这里的 "CD3E" 是 T 细胞的 Marker，你可以换成别的
FeaturePlot(nora_obj, features = "CD3E")

# 7. 画小提琴图 (看表达量高低)
VlnPlot(nora_obj, features = "CD3E")