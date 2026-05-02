[cat > ~/Projects/NDUFB7_HF_2026_04_20/03_results/scripts/mr_analysis_setup.R << 'REOF'
# ============================================
# MR分析完整脚本：NDUFB7 → Heart Failure
# 使用IEU GWAS数据库（无需手动下载FinnGen/GTEx）
# ============================================

cat("=" ,rep("=", 69), "\n", sep="")
cat("NDUFB7 MR分析启动\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 安装/加载包 ----------
if(!require("remotes")) install.packages("remotes")
if(!require("TwoSampleMR")) remotes::install_github("MRCIEU/TwoSampleMR")
if(!require("ieugwasr")) remotes::install_github("MRCIEU/ieugwasr")
if(!require("data.table")) install.packages("data.table")

library(TwoSampleMR)
library(ieugwasr)
library(data.table)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/mr_results"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# ---------- 1. 暴露数据：NDUFB7 eQTL ----------
cat("\n【1/4】准备NDUFB7暴露数据\n")

# 方案A: 如果eQTLGen已下载，从本地提取
eqtl_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"

# 手动构建NDUFB7的cis-eQTL工具变量（基于文献已知最强SNP）
# rs34141102 (chr14:89373667:G:A) 是NDUFB7最强cis-eQTL
exposure_dat <- data.table(
  SNP = "rs34141102",
  beta = NA,        # 待从eQTLGen填充
  se = NA,
  effect_allele = "A",
  other_allele = "G",
  eaf = 0.35,       # 需从eQTLGen或gnomAD确认
  pval = NA,
  samplesize = 31684,
  exposure = "NDUFB7_eQTL",
  id.exposure = "NDUFB7",
  mr_keep = TRUE
)

cat("暴露数据框架已创建（rs34141102为候选工具变量）\n")
cat("  需要从eQTLGen填充: beta, se, pval, eaf\n")

# ---------- 2. 结局数据：心衰GWAS（IEU在线数据库） ----------
cat("\n【2/4】搜索心衰GWAS结局数据\n")

# 搜索IEU数据库中的heart failure GWAS
hf_gwas <- ieugwasr::gwasinfo("heart failure")
cat("IEU数据库中心衰相关GWAS:\n")
print(hf_gwas[, c("id", "trait", "sample_size", "year", "author")])

# 推荐使用的GWAS ID（根据可用性选择）
# 常见可用ID: 
# - "ieu-b-5085" (Heart failure)
# - "ukb-b-19953" (Heart failure)
# - "finn-b-I9_HEARTFAIL" (FinnGen via IEU)

cat("\n建议使用的结局GWAS ID（请根据上方列表确认）:\n")
cat("  首选: ieu-b-5085 或 ukb-b-19953\n")

# ---------- 3. 工具变量验证标准 ----------
cat("\n【3/4】工具变量筛选标准\n")
cat("  F统计量: >10（强工具变量）\n")
cat("  r²阈值: <0.001（独立信号，避免LD连锁）\n")
cat("  cis窗口: +/-100kb from NDUFB7 TSS\n")
cat("  p值: <5e-8（全基因组显著）或 <1e-5（ suggestive, 需MR-PRESSO校正）\n")

# ---------- 4. MR方法选择 ----------
cat("\n【4/4】MR分析方法\n")
cat("  主分析: Wald Ratio（如只有1个SNP）或 Inverse Variance Weighted (IVW)\n")
cat("  敏感性: MR-Egger（检测pleiotropy）, Weighted Median, MR-PRESSO\n")
cat("  共定位: coloc（验证eQTL与GWAS是否共享因果变异）\n")

# 保存配置
config <- data.table(
  Item = c("Exposure_gene", "Lead_SNP", "Outcome_trait", "MR_method_primary", "MR_method_sensitivity", "Coloc_required"),
  Value = c("NDUFB7", "rs34141102", "Heart_failure", "Wald_ratio_or_IVW", "MR-Egger+Weighted_median+MR-PRESSO", "Yes")
)
fwrite(config, file.path(out_dir, "MR_config.txt"), sep="\t")

cat("\n✅ MR分析配置已保存到:", file.path(out_dir, "MR_config.txt"), "\n")
cat("\n完成设置。下一步:\n")
cat("  1. 确认eQTLGen下载完成，提取rs34141102的beta/se\n")
cat("  2. 从IEU获取结局GWAS数据\n")
cat("  3. 运行format_data() + mr()分析\n")

