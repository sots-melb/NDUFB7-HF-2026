[cat > mr_wald_ratio_gtex.R << 'EOF'
library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("GTEx Heart LV单SNP Wald Ratio\n")
cat("=" ,rep("=", 69), "\n", sep="")

# GTEx暴露
exp_beta <- -0.0803
exp_se <- 0.0166
exp_ea <- "T"
exp_oa <- "C"
snp <- "rs8103021"

cat("【GTEx暴露】\n")
cat("  SNP:", snp, "\n")
cat("  beta:", exp_beta, "(T↓NDUFB7)\n")
cat("  se:", exp_se, "\n")

# 结局数据（手动填入）
outcome_beta <- NA
outcome_se <- NA
outcome_ea <- NA

if(!is.na(outcome_beta)) {
  flip <- ifelse(outcome_ea == exp_ea, 1, -1)
  beta_MR <- (outcome_beta * flip) / exp_beta
  se_MR <- abs(outcome_se / exp_beta)
  pval_MR <- 2 * pnorm(-abs(beta_MR / se_MR))
  
  cat("\n【GTEx MR结果】\n")
  cat("  beta_MR:", round(beta_MR, 4), "\n")
  cat("  pval:", format(pval_MR, digits=3), "\n")
  
  if(pval_MR < 0.05) {
    cat(beta_MR > 0 ? "  🔴 心脏NDUFB7↓ → HF风险↑\n" : "  🟢 心脏NDUFB7↓ → HF风险↓\n")
  }
} else {
  cat("\n⏳ 等待rs8103021结局数据\n")
}

cat("\n完成\n")
