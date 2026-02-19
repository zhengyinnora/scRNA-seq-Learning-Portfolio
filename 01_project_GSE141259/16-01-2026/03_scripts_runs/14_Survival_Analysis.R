# ==============================================================================
# Script: 14_Survival_Analysis.R
# Purpose: 用人类 IPF 临床数据验证 ADI 特征的预后价值
# ==============================================================================

# 1. 准备工作
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
if (!requireNamespace("survminer", quietly = TRUE)) install.packages("survminer")

library(GEOquery)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# 2. 从 GEO 下载并指定保存路径 (存入 00_data_raw)
# ------------------------------------------------------------------------------
# 加上 destdir 参数，下次运行发现有文件就不会重新下载了！
message("🌐 正在抓取 GSE28042 (保存在 00_data_raw)...")
gse <- getGEO("GSE28042", destdir = "lianxi/00_data_raw", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse_data <- gse[[1]] 

# 3. 清洗临床表型数据 (提取生死簿)
# ------------------------------------------------------------------------------
message("📊 正在提取临床生死信息...")
pheno <- pData(gse_data)
clinical_data <- pheno %>%
  select(geo_accession, characteristics_ch1.1, characteristics_ch1.2) %>%
  mutate(
    Event = ifelse(grepl("Dead", characteristics_ch1.1), 2, 1),
    Time = as.numeric(gsub(".*: ", "", characteristics_ch1.2))
  )

# 4. 提取并清洗基因表达矩阵
# ------------------------------------------------------------------------------
message("🧬 正在处理基因表达矩阵...")
exprs_matrix <- exprs(gse_data)
feature_data <- fData(gse_data)

valid_genes <- feature_data$`Gene symbol` != ""
exprs_matrix <- exprs_matrix[valid_genes, ]
rownames(exprs_matrix) <- feature_data$`Gene symbol`[valid_genes]

exprs_df <- as.data.frame(exprs_matrix)
exprs_df$Gene <- rownames(exprs_matrix)
exprs_clean <- exprs_df %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean)) %>%
  tibble::column_to_rownames("Gene")

# ==============================================================================
# 🌟 下面是见证奇迹的时刻：算分与生存分析
# ==============================================================================

message("🧮 正在计算病人的 Krt8 ADI 风险得分...")

# 5. 定义我们的“恶人名单” (注意：人类的基因名全是大写！)
# 这是你前几天辛辛苦苦从老鼠里挖出来的核心机制基因
adi_signature <- c("KRT8", "SOX9", "ATF6", "TWIST1", "SPP1", "ICAM1")

# 找出存在于这个人类测序数据里的基因
available_genes <- adi_signature[adi_signature %in% rownames(exprs_clean)]

# 6. 给每个病人打分 (计算这几个基因的平均表达量)
sig_exprs <- exprs_clean[available_genes, ]
patient_scores <- colMeans(sig_exprs, na.rm = TRUE)

# 把分数加回临床数据表
clinical_data$ADI_Score <- patient_scores[clinical_data$geo_accession]

# 7. 划分高低风险组 (以全体病人的中位数为界线)
median_score <- median(clinical_data$ADI_Score, na.rm = TRUE)
clinical_data$Risk_Group <- ifelse(clinical_data$ADI_Score >= median_score, "High_ADI_Signature", "Low_ADI_Signature")

# 8. 绘制 Kaplan-Meier 生存曲线
message("🎨 正在绘制生存曲线...")
fit <- survfit(Surv(Time, Event) ~ Risk_Group, data = clinical_data)

p_surv <- ggsurvplot(fit,
                     data = clinical_data,
                     pval = TRUE,              # 显示 P 值 (小于0.05就算发财了)
                     risk.table = TRUE,        # 底部显示还在活着的病人数量
                     palette = c("red", "blue"), # 红色高风险，蓝色低风险
                     title = "IPF Survival based on Mouse ADI Signature",
                     xlab = "Time (Months)",
                     ylab = "Survival Probability",
                     legend.title = "Patient Group")

# 保存图片 (注意这里保存的是 p_surv$plot，因为 ggsurvplot 返回的是个复杂对象)
ggsave("lianxi/04_output_plots/Survival_KM_Plot.png", p_surv$plot, width = 7, height = 7)

message("🎉 大功告成！快去 04_output_plots 看看那张图！重点看图上的 P 值！")