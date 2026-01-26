# ==============================================================================
# 重新运行：第一步 (静音版)
# ==============================================================================

# 1. 安装构建工具
install.packages("remotes")

# 2. 安装 BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 3. 安装依赖包 (binary)
difficult_pkgs <- c("sf", "terra", "leidenbase", "igraph")
install.packages(difficult_pkgs, type = "binary")

# 4. 【关键修改】安装 Bioconductor 依赖 (强制不更新，不问问题)
BiocManager::install(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats",
                       "limma", "lme4", "S4Vectors", "SingleCellExperiment",
                       "SummarizedExperiment", "batchelor", "HDF5Array",
                       "terra", "ggrastr"),
                     update = FALSE, ask = FALSE) # <--- 重点在这里

message("🎉 第一步：地基铺设完毕！")

# ==============================================================================
# Part 2: Installing Monocle3 (静音版 - Silent Mode)
# ==============================================================================

message("开始从 GitHub 下载 Monocle3... (这取决于网速，可能需要几分钟 🐢->🐇)")

# 1. 强制安装已下架的 grr 包
remotes::install_version("grr", version = "0.9.5", repos = "https://cloud.r-project.org")
# upgrade = "never": 
# 意思是：“只管装 Monocle3，别管我电脑里其他包旧不旧，千万别问我要不要更新！直接装！”
remotes::install_github('cole-trapnell-lab/monocle3', upgrade = "never")

# 验证环节：尝试加载
library(monocle3)
message("🎉 恭喜！Monocle3 安装并加载成功！任务完成！")