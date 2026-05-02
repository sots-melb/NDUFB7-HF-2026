library(data.table)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/mr_results"
cat("=" ,rep("=", 69), "\n", sep="")
cat("MR Wald Ratio: 单SNP快速因果估计\n")
cat("=" ,rep("=", 69), "\n", sep="")

# 读取最强SNP
exp <- fread(file.path(out_dir, "NDUFB7_exposure_eQTLGen_strongIV.csv"))
top <- exp[which.min(pval)]
cat("最强工具变量:\n")
print(top[, c("SNP","beta","se","pval","effect_allele","other_allele")])

# Wald Ratio公式: beta_outcome / beta_exposure
# 需要结局GWAS中该SNP的beta/se
# 由于gwasinfo API故障，使用手动已知值或IEU在线查询

cat("\n【手动Wald Ratio计算】\n")
cat("公式: b_MR = b_outcome / b_exposure\n")
cat("       se_MR = se_outcome / |b_exposure|\n")

# 假设从IEU数据库查到的rs11085898在心衰GWAS中的效应（示例值，需替换）
# 实际应从 outcome GWAS 提取，这里展示计算框架
b_exp <- top$beta
se_exp <- top$se
p_exp <- top$pval

cat("\n暴露 (eQTLGen):\n")
cat("  rs11085898: beta=", round(b_exp,4), ", se=", round(se_exp,4), ", p=", format(p_exp, digits=2), "\n")

cat("\n⚠️ 需要结局GWAS中rs11085898的beta/se才能计算Wald Ratio\n")
cat("   方案A: 等待IEU API恢复，用extract_outcome_data()\n")
cat("   方案B: 手动查询 https://gwas.mrcieu.ac.uk/ 获取rs11085898在心衰GWAS中的值\n")
cat("   方案C: 使用已发表的FinnGen/HERMES汇总统计本地提取\n")

# 保存框架供后续填入
template <- data.table(
  SNP = "rs11085898",
  beta_exposure = b_exp,
  se_exposure = se_exp,
  pval_exposure = p_exp,
  beta_outcome = NA,
  se_outcome = NA,
  pval_outcome = NA,
  beta_MR = NA,
  se_MR = NA,
  pval_MR = NA
)
fwrite(template, file.path(out_dir, "Wald_ratio_template.txt"), sep="\t")

cat("\n✅ Wald Ratio模板已保存，填入结局数据后自动计算\n")
cat("   文件:", file.path(out_dir, "Wald_ratio_template.txt"), "\n")

cat("\n完成\n")
