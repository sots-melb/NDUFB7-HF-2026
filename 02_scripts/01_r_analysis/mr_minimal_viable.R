library(data.table)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/mr_results"
cat("=" ,rep("=", 69), "\n", sep="")
cat("MR最小可行方案: eQTLGen暴露 + 手动结局查询\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 1. 暴露数据（已完成） ----------
exp <- fread(file.path(out_dir, "NDUFB7_exposure_eQTLGen_strongIV.csv"))
top <- exp[which.min(pval)]
cat("最强工具变量:\n")
print(top[, c("SNP","SNPChr","SNPPos","beta","se","pval","effect_allele","other_allele")])

# ---------- 2. 手动查询结局GWAS（指导用户） ----------
cat("\n【手动查询指南】\n")
cat("你需要在本地浏览器查询以下SNP在心衰GWAS中的效应:\n\n")

snps_to_query <- head(exp[order(pval)]$SNP, 5)
for(snp in snps_to_query) {
  cat("  ", snp, "\n")
}

cat("\n查询方法:\n")
cat("  方法1: IEU GWAS数据库\n")
cat("    访问: https://gwas.mrcieu.ac.uk/region.html\n")
cat("    输入SNP，选择结局GWAS（如ukb-b-19953或finn-b-I9_HEARTFAIL）\n")
cat("    记录: beta_outcome, se_outcome, pval_outcome, effect_allele, other_allele\n\n")

cat("  方法2: PhenoScanner\n")
cat("    访问: http://www.phenoscanner.medschl.cam.ac.uk/\n")
cat("    输入SNP，选择Heart failure相关traits\n\n")

cat("  方法3: GWAS Catalog\n")
cat("    访问: https://www.ebi.ac.uk/gwas/\n")
cat("    搜索SNP，查看关联studies\n\n")

# ---------- 3. 本地Wald Ratio计算框架 ----------
cat("\n【Wald Ratio计算框架】\n")
b_exp <- top$beta
se_exp <- top$se

cat("暴露效应 (eQTLGen rs11085898):\n")
cat("  beta_exp =", round(b_exp, 4), "\n")
cat("  se_exp   =", round(se_exp, 4), "\n")

cat("\n填入结局效应后计算:\n")
cat("  beta_MR  = beta_outcome / beta_exp\n")
cat("  se_MR    = se_outcome / |beta_exp|\n")
cat("  p_MR     = 2 * pnorm(-abs(beta_MR / se_MR))\n")
cat("  OR       = exp(beta_MR)\n")
cat("  CI_lower = exp(beta_MR - 1.96 * se_MR)\n")
cat("  CI_upper = exp(beta_MR + 1.96 * se_MR)\n")

# 保存查询模板
template <- data.table(
  SNP = snps_to_query,
  beta_exp = exp[match(snps_to_query, SNP)]$beta,
  se_exp = exp[match(snps_to_query, SNP)]$se,
  beta_outcome = NA,
  se_outcome = NA,
  pval_outcome = NA,
  outcome_source = NA,
  beta_MR = NA,
  se_MR = NA,
  pval_MR = NA,
  OR = NA,
  CI_lower = NA,
  CI_upper = NA
)
fwrite(template, file.path(out_dir, "MR_manual_query_template.txt"), sep="\t")

cat("\n✅ 查询模板已保存:", file.path(out_dir, "MR_manual_query_template.txt"), "\n")
cat("   填入结局数据后，用Excel或R计算Wald Ratio\n")

cat("\n完成\n")
