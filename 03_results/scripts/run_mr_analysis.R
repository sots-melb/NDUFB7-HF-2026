
library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("TwoSampleMR: NDUFB7 -> Heart Failure 分析\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 安装 ----------
if(!require("remotes")) install.packages("remotes", repos="https://cloud.r-project.org")
if(!require("TwoSampleMR")) {
  cat("安装TwoSampleMR...\n")
  remotes::install_github("MRCIEU/TwoSampleMR", force=FALSE)
}
if(!require("ieugwasr")) {
  cat("安装ieugwasr...\n")
  remotes::install_github("MRCIEU/ieugwasr", force=FALSE)
}

library(TwoSampleMR)
library(ieugwasr)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/mr_results"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# ---------- 1. 读取暴露数据 ----------
cat("\n【1/5】读取NDUFB7暴露数据\n")
exp_file <- file.path(out_dir, "NDUFB7_exposure_eQTLGen_strongIV.csv")
if(!file.exists(exp_file)) {
  exp_file <- file.path(out_dir, "NDUFB7_exposure_eQTLGen_all38.csv")
}

exposure_raw <- fread(exp_file)
cat("暴露SNP数:", nrow(exposure_raw), "\n")

# 格式化为TwoSampleMR格式
exposure_dat <- format_data(
  exposure_raw,
  type = "exposure",
  snps = NULL,
  header = TRUE,
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  samplesize_col = "samplesize",
  phenotype_col = "exposure",
  id_col = "id.exposure"
)

cat("格式化后暴露数据:\n")
print(head(exposure_dat[, c("SNP","beta","se","pval","effect_allele","other_allele")]))

# 筛选p<5e-8（严格）或p<5e-6（宽松，需MR-PRESSO）
exposure_dat <- clump_data(exposure_dat, clump_p1=5e-6, clump_p2=5e-6, clump_r2=0.001, clump_kb=10000)
cat("Clumping后SNP数:", nrow(exposure_dat), "\n")

# ---------- 2. 搜索结局GWAS ----------
cat("\n【2/5】搜索IEU GWAS数据库: Heart Failure\n")
hf_search <- ieugwasr::gwasinfo("heart failure")
cat("找到", nrow(hf_search), "个相关GWAS\n")

# 显示前10个
if(nrow(hf_search) > 0) {
  print(hf_search[1:min(10,nrow(hf_search)), c("id","trait","sample_size","year","author","population")])
  
  # 选择最佳结局
  # 优先: UKB or FinnGen via IEU, 欧洲人群, 大样本
  preferred <- hf_search[grepl("UKB|FinnGen|HERMES|heart failure", hf_search$trait, ignore.case=TRUE), ]
  cat("\n优先候选:\n")
  print(preferred[, c("id","trait","sample_size","year")])
  
  # 保存候选列表
  fwrite(hf_search, file.path(out_dir, "HF_GWAS_candidates.txt"), sep="\t")
}

# ---------- 3. 获取结局数据（手动指定ID） ----------
cat("\n【3/5】获取结局GWAS数据\n")

# 尝试多个候选ID
candidate_ids <- c("ieu-b-5085", "ukb-b-19953", "finn-b-I9_HEARTFAIL", "ebi-a-GCST009119")

outcome_dat <- NULL
for(oid in candidate_ids) {
  cat("尝试:", oid, "... ")
  tryCatch({
    od <- extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = oid,
      proxies = TRUE,
      rsq = 0.8
    )
    if(!is.null(od) && nrow(od) > 0) {
      outcome_dat <- od
      cat("✅ 成功 (", nrow(od), " SNP匹配)\n")
      break
    } else {
      cat("⚠️ 无匹配SNP\n")
    }
  }, error = function(e) {
    cat("❌ 失败:", conditionMessage(e), "\n")
  })
}

if(is.null(outcome_dat) || nrow(outcome_dat) == 0) {
  cat("\n❌ 所有候选ID均失败。保存候选列表供手动选择。\n")
  cat("建议操作: 访问 https://gwas.mrcieu.ac.uk/ 搜索 'heart failure'\n")
  cat("或查看保存的 HF_GWAS_candidates.txt\n")
  
  # 保存暴露数据供后续使用
  fwrite(exposure_dat, file.path(out_dir, "exposure_formatted.txt"), sep="\t")
  cat("暴露数据已保存，可后续手动匹配结局\n")
} else {
  # ---------- 4.  harmonize ----------
  cat("\n【4/5】Harmonize暴露与结局\n")
  dat <- harmonise_data(exposure_dat, outcome_dat, action=2)
  cat("Harmonize后:", nrow(dat), "SNP\n")
  
  # ---------- 5. MR分析 ----------
  cat("\n【5/5】MR分析\n")
  
  # 主分析
  res <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  cat("\nMR结果:\n")
  print(res[, c("method","nsnp","b","se","pval")])
  
  # 敏感性
  cat("\n异质性检验:\n")
  het <- mr_heterogeneity(dat)
  print(het[, c("method","Q","Q_df","Q_pval")])
  
  cat("\n多效性检验 (MR-Egger intercept):\n")
  pleio <- mr_pleiotropy_test(dat)
  print(pleio[, c("egger_intercept","se","pval")])
  
  # 保存
  fwrite(res, file.path(out_dir, "MR_results.txt"), sep="\t")
  fwrite(dat, file.path(out_dir, "MR_harmonized_data.txt"), sep="\t")
  cat("\n✅ MR结果已保存到:", out_dir, "\n")
}

cat("\n完成\n")
