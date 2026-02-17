# ==============================================================================
# Script: 12_Install_NicheNet.R
# Purpose: 安装 NicheNet 包 (devtools 方式)
# ==============================================================================

# 1. 安装 devtools (如果还没有)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# 2. 安装 nichenetr (从 GitHub)
#    这是一个纯 R 包，通常安装比较顺利
message("⬇️ 正在安装 nichenetr...")

if (!requireNamespace("nichenetr", quietly = TRUE)) {
  tryCatch({
    # 尝试镜像安装
    devtools::install_git("https://gitclone.com/github.com/saeyslab/nichenetr.git")
    message("✅ nichenetr 安装成功！")
  }, error = function(e) {
    message("⚠️ 镜像失败，尝试官方源...")
    devtools::install_github("saeyslab/nichenetr")
  })
}

# 3. 安装 tidyverse (数据处理神器，NicheNet 极其依赖它)
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

message("🎉 软件包准备就绪！请确保你已经手动下载了 3 个 .rds 文件并上传到了 lianxi 文件夹！")