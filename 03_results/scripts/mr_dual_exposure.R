library(data.table)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/mr_results"
cat("=" ,rep("=", 69), "\n", sep="")
cat("MR双暴露策略: eQTLGen(全血) + GTEx(心脏)\n")
cat("=" ,rep("=", 69), "\n", sep="")

# 暴露1: eQTLGen全血 (38 SNP)
exp_blood <- fread(file.path(out_dir, "NDUFB7_exposure_eQTLGen_strongIV.csv"))
cat("【暴露1】eQTLGen全血:", nrow(exp_blood), "SNP\n")

# 暴露2: GTEx Heart LV (1 SNP, 从eGene提取)
exp_heart <- data.table(
  SNP = "rs8103021",
  beta = -0.0803,      # GTEx slope (标准化)
  se = 0.0166,
  pval = 6.55e-6,
  effect_allele = "T",
  other_allele = "C",
  samplesize = 387,    # GTEx Heart LV样本量
  exposure = "NDUFB7_eQTL_HeartLV",
  id.exposure = "NDUFB7_Heart",
  Gene = "NDUFB7",
  GeneChr = 19,
  GenePos = 14641837
)
cat("【暴露2】GTEx Heart LV: 1 SNP (rs8103021)\n")

cat("\n【MR策略选择】\n")
cat("  方案A: 主分析用eQTLGen(38 SNP, IVW) — 统计效力最强\n")
cat("  方案B: 敏感性分析用GTEx(1 SNP, Wald Ratio) — 组织特异性验证\n")
cat("  方案C: 若两暴露方向一致，增强因果推断; 若不一致，提示组织依赖性\n")

cat("\n【当前阻塞】\n")
cat("  需手动查询结局GWAS中rs11085898和rs8103021的效应\n")
cat("  或等待IEU API恢复自动提取\n")

# 保存双暴露配置
fwrite(exp_heart, file.path(out_dir, "NDUFB7_exposure_GTEx_HeartLV.txt"), sep="\t")
cat("\n✅ GTEx暴露已保存\n")

cat("\n【下一步】\n")
cat("  1. 手动查询PhenoScanner/GWAS Catalog获取结局效应\n")
cat("  2. 或尝试IEU extract_outcome_data()\n")
cat("  3. 计算Wald Ratio: beta_outcome / beta_exposure\n")

cat("\n完成\n")
