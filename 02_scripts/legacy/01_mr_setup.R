#!/usr/bin/env Rscript
setwd("~/Projects/NDUFB7_HF_{2026_04_20}")

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     TwoSampleMR安装与MR数据准备                          ║\n")
cat("║     目标: NDUFB7 eQTL → HF GWAS 因果推断                 ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# --- 0. 安装TwoSampleMR ---
cat("========== 0. 安装TwoSampleMR ==========\n")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
}

install_success <- FALSE
tryCatch({
  remotes::install_github("MRCIEU/TwoSampleMR", force = FALSE, quiet = TRUE)
  install_success <- TRUE
  cat("✅ TwoSampleMR安装成功\n")
}, error = function(e) {
  cat("❌ TwoSampleMR安装失败:", conditionMessage(e), "\n")
  cat("备选: 使用本地MR基础方法（wald ratio, IVW公式手动计算）\n")
})

if (install_success) {
  library(TwoSampleMR)
  cat("✅ TwoSampleMR版本:", packageVersion("TwoSampleMR"), "\n")
}

# --- 1. 数据下载策略 ---
cat("\n========== 1. MR数据下载策略 ==========\n")

mr_plan <- data.frame(
  数据 = c("eQTLGen血eQTL", "GTEx心脏eQTL", "FinnGen心衰GWAS", "UKBB心衰GWAS"),
  用途 = c("暴露(NDUFB7)", "暴露验证(心脏特异性)", "结局(HF)", "结局验证"),
  状态 = c("待下载", "待下载", "待下载", "备选"),
  下载URL = c(
    "https://www.eqtlgen.org/phase1.html",
    "https://gtexportal.org/home/datasets",
    "https://www.finngen.fi/en/access_results",
    "https://www.nealelab.is/uk-biobank"
  ),
  文件格式 = c("txt.gz", "tsv.gz", "tsv.gz", "tsv.gz"),
  备注 = c(
    "需注册或公开下载，搜索NDUFB7 rsID",
    "GTEx v8 Heart - Left Ventricle eQTL",
    "FinnGen R10 Heart failure",
    "需申请，数据量大"
  )
)

print(mr_plan, row.names = FALSE)

# --- 2. 生成下载命令 ---
cat("\n========== 2. 生成下载命令 ==========\n")

download_cmds <- '
#!/bin/bash
# MR数据下载脚本（手动执行或放入后台）
# 建议逐个下载，避免网络拥堵

mkdir -p ~/Projects/NDUFB7_HF_{2026_04_20}/01_data/03_mr_data

# 1. eQTLGen (需从网站手动下载，无直接FTP)
# 访问: https://www.eqtlgen.org/phase1.html
# 下载: cis-eQTL summary statistics
# 搜索基因: NDUFB7 (chr19: 需确认具体位置)

# 2. GTEx心脏eQTL (可直接wget)
cd ~/Projects/NDUFB7_HF_{2026_04_20}/01_data/03_mr_data
wget -c https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_all_associations/Heart_Left_Ventricle.v8.allpairs.txt.gz
echo "GTEx心脏eQTL下载完成"

# 3. FinnGen心衰GWAS (需从网站下载，或申请API)
# 访问: https://www.finngen.fi/en/access_results
# 下载: Endpoint: I9_HEARTFAIL (Heart failure)

# 4. 备选: UK Biobank GWAS summary
# 需通过GWAS Catalog或Neale Lab获取

echo "MR数据下载脚本准备完成"
'

cat(download_cmds)
writeLines(download_cmds, "02_scripts/06_mr_preparation/02_download_mr.sh")
Sys.chmod("02_scripts/06_mr_preparation/02_download_mr.sh", mode = "0755")
cat("\n✅ 下载脚本已保存: 02_scripts/06_mr_preparation/02_download_mr.sh\n")

# --- 3. 本地MR分析框架（即使TwoSampleMR安装失败也可用） ---
cat("\n========== 3. 本地MR基础分析框架 ==========\n")

local_mr_code <- '
# 本地MR分析（ Wald Ratio + IVW 手动实现）
# 适用于已获取的eQTL和GWAS summary数据

library(dplyr)

# 读取eQTL数据（NDUFB7的cis-eQTL）
eqtl <- read.table("eqtl_file.txt.gz", header = TRUE, sep = "\t")

# 筛选NDUFB7的显著SNP（p < 5e-8）
ndufb7_snps <- eqtl %>% 
  filter(gene == "NDUFB7", pval < 5e-8) %>%
  select(SNP, beta, se, pval, effect_allele, other_allele)

# 读取GWAS数据（HF）
gwas <- read.table("gwas_file.txt.gz", header = TRUE, sep = "\t")

# 合并
merged <- inner_join(ndufb7_snps, gwas, by = "SNP")

# Wald Ratio for each SNP
merged$wr_beta <- merged$beta.y / merged$beta.x
merged$wr_se <- merged$se.y / abs(merged$beta.x)

# IVW meta-analysis
ivw_beta <- sum(merged$wr_beta / merged$wr_se^2) / sum(1 / merged$wr_se^2)
ivw_se <- sqrt(1 / sum(1 / merged$wr_se^2))
ivw_p <- 2 * pnorm(-abs(ivw_beta / ivw_se))

cat("IVW结果: beta =", ivw_beta, "se =", ivw_se, "p =", ivw_p, "\\n")
'

cat(local_mr_code)
writeLines(local_mr_code, "02_scripts/06_mr_preparation/03_local_mr_analysis.R")
cat("✅ 本地MR框架已保存: 02_scripts/06_mr_preparation/03_local_mr_analysis.R\n")

# --- 4. 记录MR计划 ---
cat("\n========== 4. MR执行计划 ==========\n")
cat("【Phase 1】数据获取（本周）\n")
cat("  □ eQTLGen: 访问网站，下载cis-eQTL，提取NDUFB7的top SNP\n")
cat("  □ GTEx: 执行wget下载Heart_LV eQTL\n")
cat("  □ FinnGen: 申请或下载HF GWAS summary\n")
cat("\n【Phase 2】分析执行（下周）\n")
cat("  □ 数据清洗: 统一SNP命名，匹配allele，计算palindromic SNP\n")
cat("  □ MR分析: Wald Ratio, IVW, MR-Egger, Weighted Median\n")
cat("  □ 敏感性: 留一法, Cochran Q, MR-PRESSO, 漏斗图\n")
cat("  □ 结果: 因果效应估计, 方向性验证, 发表级图表\n")

# 保存计划
write.csv(mr_plan, "03_results/18_mr_data/60_mr_data_plan.csv", row.names = FALSE)
cat("\n✅ MR准备完成，计划保存至 03_results/18_mr_data/\n")
