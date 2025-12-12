install.packages("Seurat")
install.packages("tidyverse")
library(Seurat)
library(tidyverse)
# 设置CRAN镜像（推荐用德国镜像，快又稳）
options(repos = c(CRAN = "https://cloud.r-project.org"))
# 安装Seurat和常用可视化包
install.packages(c("Seurat", "patchwork", "ggplot2", "dplyr"))
install.packages(c("Seurat", "patchwork", "ggplot2", "dplyr"))
nora_obj <- pbmc_small
nora_obj
DimPlot(nora_obj)
FeaturePlot(nora_obj, features = "CD3E")
VlnPlot(nora_obj, features = "CD3E")
# 看一下对象结构nora_obj
# 看有哪些基因
rownames(nora_obj)
# 按照cluster画FeaturePlot
FeaturePlot(nora_obj, features = "CD3E")   # 看T细胞marker
# 看一下对象结构
nora_obj
# 看有哪些基因
rownames(nora_obj)
# 按照cluster画FeaturePlot
FeaturePlot(nora_obj, features = "CD3E")   # 看T细胞marker
# 看一下对象结构
nora_obj
# 看有哪些基因
rownames(nora_obj)
# 按照cluster画FeaturePlot
FeaturePlot(nora_obj, features = "CD3E")   # 看T细胞marker
