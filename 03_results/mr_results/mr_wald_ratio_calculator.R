[cd ~/Projects/NDUFB7_HF_2026_04_20/03_results/mr_results

cat > mr_wald_ratio_calculator.R << 'EOF'
library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("MR Wald Ratio计算器\n")
cat("填入手动查询的结局数据后运行\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ========== 暴露数据（已就绪） ==========
exp <- fread("NDUFB7_exposure_eQTLGen_strongIV.csv")
top <- exp[which.min(pval)]

cat("\n【暴露】\n")
cat("  SNP:", top$SNP, "\n")
cat("  暴露beta:", round(top$beta, 4), "\n")
cat("  暴露se:", round(top$se, 4), "\n")
cat("  暴露pval:", format(top$pval, digits=2), "\n")
cat("  效应等位基因:", top$effect_allele, "\n")

# ========== 结局数据（手动填入） ==========
# 请根据PhenoScanner/GWAS Catalog查询结果填入以下数值
# 注意：如果结局的effect_allele与暴露不同，需要flip beta_outcome

outcome_beta <- NA    # 填入结局beta（或ln(OR)）
outcome_se <- NA      # 填入结局se
outcome_pval <- NA    # 填入结局pval
outcome_ea <- NA      # 填入结局effect_allele
outcome_oa <- NA      # 填入结局other_allele
outcome_n <- NA       # 填入结局样本量

cat("\n【结局】（当前为NA，请手动填入）\n")
cat("  结局beta:", outcome_beta, "\n")
cat("  结局se:", outcome_se, "\n")
cat("  结局pval:", outcome_pval, "\n")

# ========== 自动计算 ==========
if(!is.na(outcome_beta) && !is.na(outcome_se)) {
  
  # 检查等位基因方向
  flip <- 1
  if(!is.na(outcome_ea) && outcome_ea != top$effect_allele) {
    if(outcome_ea == top$other_allele) {
      flip <- -1
      cat("\n⚠️ 结局effect_allele与暴露相反，自动flip\n")
    } else {
      cat("\n❌ 等位基因不匹配，无法计算\n")
      quit(status=1)
    }
  }
  
  # Wald Ratio
  beta_MR <- (outcome_beta * flip) / top$beta
  se_MR <- abs(outcome_se / top$beta)
  pval_MR <- 2 * pnorm(-abs(beta_MR / se_MR))
  OR <- exp(beta_MR)
  CI_lower <- exp(beta_MR - 1.96 * se_MR)
  CI_upper <- exp(beta_MR + 1.96 * se_MR)
  
  cat("\n【MR结果】\n")
  cat("  Wald Ratio beta:", round(beta_MR, 4), "\n")
  cat("  SE:", round(se_MR, 4), "\n")
  cat("  P-value:", format(pval_MR, digits=3), "\n")
  cat("  OR:", round(OR, 3), "\n")
  cat("  95% CI:", round(CI_lower, 3), "-", round(CI_upper, 3), "\n")
  
  # 解释
  cat("\n【解读】\n")
  if(pval_MR < 0.05) {
    if(beta_MR > 0) {
      cat("  🔴 遗传预测的NDUFB7表达增加与心衰风险增加相关\n")
    } else {
      cat("  🟢 遗传预测的NDUFB7表达增加与心衰风险降低相关（保护性）\n")
    }
  } else {
    cat("  ⚪ 未发现NDUFB7表达与心衰的显著因果关联 (p>0.05)\n")
  }
  
  # 保存
  result <- data.table(
    SNP = top$SNP,
    exposure = "NDUFB7_eQTL",
    outcome = "Heart_failure",
    beta_MR = beta_MR,
    se_MR = se_MR,
    pval_MR = pval_MR,
    OR = OR,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    method = "Wald_Ratio"
  )
  fwrite(result, "MR_wald_ratio_result.txt", sep="\t")
  cat("\n✅ 结果已保存: MR_wald_ratio_result.txt\n")
  
} else {
  cat("\n⏳ 等待结局数据填入...\n")
  cat("  请编辑本脚本，填入outcome_beta, outcome_se, outcome_pval\n")
  cat("  然后重新运行: Rscript mr_wald_ratio_calculator.R\n")
}

cat("\n完成\n")
