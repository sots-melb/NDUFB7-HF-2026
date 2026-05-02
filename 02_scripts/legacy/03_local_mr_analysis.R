
# 本地MR分析（ Wald Ratio + IVW 手动实现）
# 适用于已获取的eQTL和GWAS summary数据

library(dplyr)

# 读取eQTL数据（NDUFB7的cis-eQTL）
eqtl <- read.table("eqtl_file.txt.gz", header = TRUE, sep = "	")

# 筛选NDUFB7的显著SNP（p < 5e-8）
ndufb7_snps <- eqtl %>% 
  filter(gene == "NDUFB7", pval < 5e-8) %>%
  select(SNP, beta, se, pval, effect_allele, other_allele)

# 读取GWAS数据（HF）
gwas <- read.table("gwas_file.txt.gz", header = TRUE, sep = "	")

# 合并
merged <- inner_join(ndufb7_snps, gwas, by = "SNP")

# Wald Ratio for each SNP
merged$wr_beta <- merged$beta.y / merged$beta.x
merged$wr_se <- merged$se.y / abs(merged$beta.x)

# IVW meta-analysis
ivw_beta <- sum(merged$wr_beta / merged$wr_se^2) / sum(1 / merged$wr_se^2)
ivw_se <- sqrt(1 / sum(1 / merged$wr_se^2))
ivw_p <- 2 * pnorm(-abs(ivw_beta / ivw_se))

cat("IVW结果: beta =", ivw_beta, "se =", ivw_se, "p =", ivw_p, "\n")

