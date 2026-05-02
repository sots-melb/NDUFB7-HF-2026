
library(data.table)
library(TwoSampleMR)
library(ieugwasr)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/mr_results"

cat("=" ,rep("=", 69), "\n", sep="")
cat("TwoSampleMR修复版: NDUFB7 -> Heart Failure\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 1. 手动构建暴露数据（绕过format_data） ----------
cat("\n【1/5】手动构建暴露数据\n")
exp_raw <- fread(file.path(out_dir, "NDUFB7_exposure_eQTLGen_strongIV.csv"))
cat("读取", nrow(exp_raw), "SNP\n")

# 手动创建TwoSampleMR格式
exposure_dat <- data.frame(
  SNP = exp_raw$SNP,
  beta = as.numeric(exp_raw$beta),
  se = as.numeric(exp_raw$se),
  effect_allele = exp_raw$effect_allele,
  other_allele = exp_raw$other_allele,
  pval = as.numeric(exp_raw$pval),
  samplesize = as.numeric(exp_raw$samplesize),
  exposure = "NDUFB7_eQTL",
  id.exposure = "NDUFB7",
  eaf = 0.5,  # 假设，eQTLGen未提供eaf
  mr_keep = TRUE,
  stringsAsFactors = FALSE
)

cat("暴露数据预览:\n")
print(head(exposure_dat[, c("SNP","beta","se","pval","effect_allele","other_allele")]))

# Clumping（需要LD参考，在线使用EUR）
cat("\n执行clumping (EUR参考, r2<0.001, kb=10000)...\n")
tryCatch({
  exposure_dat <- clump_data(exposure_dat, clump_p1=5e-6, clump_p2=5e-6, clump_r2=0.001, clump_kb=10000)
  cat("Clumping后:", nrow(exposure_dat), "SNP\n")
}, error = function(e) {
  cat("⚠️ Clumping失败（可能需要LD token或本地LD）:", conditionMessage(e), "\n")
  cat("   继续用全部SNP（假设已独立）\n")
})

# ---------- 2. 搜索结局GWAS ----------
cat("\n【2/5】搜索IEU GWAS: Heart Failure\n")
hf_search <- ieugwasr::gwasinfo("heart failure")
cat("找到", nrow(hf_search), "个GWAS\n")

if(nrow(hf_search) > 0) {
  # 显示关键信息
  disp_cols <- intersect(c("id","trait","sample_size","year","author","population","ncase","ncontrol"), names(hf_search))
  print(hf_search[1:min(15,nrow(hf_search)), ..disp_cols])
  
  # 保存
  fwrite(hf_search, file.path(out_dir, "HF_GWAS_candidates.txt"), sep="\t")
  
  # 自动选择最佳
  # 优先：欧洲人群、样本量大、有case/control数
  best_idx <- which.max(hf_search$sample_size * grepl("European|Finnish|UK", hf_search$population, ignore.case=TRUE))
  if(length(best_idx) == 0) best_idx <- 1
  best_id <- hf_search$id[best_idx]
  cat("\n自动选择最佳结局ID:", best_id, "\n")
  cat("  Trait:", hf_search$trait[best_idx], "\n")
  cat("  Sample:", hf_search$sample_size[best_idx], "\n")
} else {
  cat("❌ 未找到心衰GWAS\n")
  best_id <- NULL
}

# ---------- 3. 获取结局数据 ----------
outcome_dat <- NULL
if(!is.null(best_id)) {
  cat("\n【3/5】获取结局数据 (ID:", best_id, ")\n")
  tryCatch({
    outcome_dat <- extract_outcome_data(
      snps = exposure_dat$SNP,
      outcomes = best_id,
      proxies = TRUE,
      rsq = 0.8
    )
    cat("结局数据:", nrow(outcome_dat), "SNP匹配\n")
    if(nrow(outcome_dat) > 0) {
      print(head(outcome_dat[, c("SNP","beta","se","pval","effect_allele","other_allele","outcome")]))
    }
  }, error = function(e) {
    cat("❌ 获取失败:", conditionMessage(e), "\n")
  })
}

# ---------- 4. Harmonize + MR ----------
if(!is.null(outcome_dat) && nrow(outcome_dat) > 0) {
  cat("\n【4/5】Harmonize\n")
  dat <- harmonise_data(exposure_dat, outcome_dat, action=2)
  cat("Harmonize后:", nrow(dat), "SNP\n")
  
  cat("\n【5/5】MR分析\n")
  
  # 选择方法
  methods <- c("mr_ivw")
  if(nrow(dat) == 1) methods <- c("mr_wald_ratio")
  if(nrow(dat) >= 3) methods <- c(methods, "mr_egger_regression", "mr_weighted_median")
  
  res <- mr(dat, method_list=methods)
  cat("\nMR结果:\n")
  print(res[, c("method","nsnp","b","se","pval")])
  
  # 敏感性
  if(nrow(dat) >= 3) {
    cat("\n异质性:\n")
    het <- mr_heterogeneity(dat)
    print(het[, c("method","Q","Q_df","Q_pval")])
    
    cat("\n多效性 (MR-Egger intercept):\n")
    pleio <- mr_pleiotropy_test(dat)
    print(pleio[, c("egger_intercept","se","pval")])
    
    cat("\nLeave-one-out:\n")
    loo <- mr_leaveoneout(dat)
    print(head(loo[, c("SNP","b","se","pval")]))
  }
  
  # 保存
  fwrite(as.data.table(res), file.path(out_dir, "MR_results.txt"), sep="\t")
  fwrite(as.data.table(dat), file.path(out_dir, "MR_harmonized_data.txt"), sep="\t")
  
  # 生成森林图数据
  cat("\n✅ MR分析完成，结果保存到:", out_dir, "\n")
  
  # 关键结论
  ivw_p <- res$pval[res$method=="Inverse variance weighted"]
  if(length(ivw_p) > 0 && !is.na(ivw_p)) {
    if(ivw_p < 0.05) {
      cat("\n🔥 关键发现: NDUFB7 cis-eQTL与心衰显著关联 (IVW p=", format(ivw_p, digits=3), ")\n")
    } else {
      cat("\n⚠️ NDUFB7 cis-eQTL与心衰无显著关联 (IVW p=", format(ivw_p, digits=3), ")\n")
    }
  }
} else {
  cat("\n❌ 无结局数据，保存暴露数据供后续手动匹配\n")
  fwrite(as.data.table(exposure_dat), file.path(out_dir, "exposure_only.txt"), sep="\t")
  cat("建议: 查看 HF_GWAS_candidates.txt 手动选择结局ID\n")
}

cat("\n完成\n")
