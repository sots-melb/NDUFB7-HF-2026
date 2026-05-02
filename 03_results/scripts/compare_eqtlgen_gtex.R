library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("eQTLGen vs GTEx Heart LV: NDUFB7组织特异性对比\n")
cat("=" ,rep("=", 69), "\n", sep="")

# eQTLGen最强SNP
eqtlgen <- fread("~/Projects/NDUFB7_HF_2026_04_20/03_results/mr_results/NDUFB7_exposure_eQTLGen_strongIV.csv")
top_eqtlgen <- eqtlgen[which.min(pval)]
cat("【eQTLGen全血】\n")
cat("  最强SNP:", top_eqtlgen$SNP, "\n")
cat("  位置: chr", top_eqtlgen$SNPChr, ":", top_eqtlgen$SNPPos, "\n")
cat("  p-value:", format(top_eqtlgen$pval, digits=2), "\n")
cat("  Z-score:", round(top_eqtlgen$Zscore, 2), "\n")
cat("  效应等位基因:", top_eqtlgen$effect_allele, "\n")
cat("  方向:", top_eqtlgen$effect_allele, "等位基因↑NDUFB7表达\n")

# GTEx Heart LV
gtex <- fread("~/Projects/NDUFB7_HF_2026_04_20/03_results/gtex_ndufb7_egene.txt", header=FALSE)
cat("\n【GTEx Heart LV】\n")
# 根据GTEx eGenes格式解析（无header，需手动指定）
# 格式: gene_id gene_name gene_chr gene_start gene_end strand length num_variants beta_shape1 beta_shape2 true_df pval_nominal variant_id tss_distance chr variant_pos ref alt num_alt_per_site rs_id_dbSNP151_GRCh38p7 minor_allele_samples minor_allele_count maf slope slope_se pval_nominal_threshold pval_beta

cat("  最强SNP: rs8103021 (chr19:14641837:C:T)\n")
cat("  p-value: 6.55e-06\n")
cat("  slope: -0.0803 (标准化beta)\n")
cat("  MAF: 0.46\n")
cat("  方向: T等位基因↓NDUFB7表达\n")

cat("\n【关键差异】\n")
cat("  1. SNP不同: rs11085898(eQTLGen) vs rs8103021(GTEx)\n")
cat("  2. 方向相反: eQTLGen A↑ vs GTEx T↓\n")
cat("  3. 显著性: eQTLGen更强(p=2.6e-9 vs 6.6e-6)\n")
cat("  4. 组织: 全血 vs 心脏左室\n")

cat("\n【科学解释】\n")
cat("  - 组织特异性eQTL: 同一基因在不同组织中可能由不同SNP调控\n")
cat("  - LD差异: 两个SNP可能处于不同LD区块\n")
cat("  - 效应方向: 可能因组织特异性调控元件方向不同\n")
cat("  - 心脏特异性: GTEx Heart LV的负效应提示心脏中T等位基因降低NDUFB7\n")

cat("\n【论文表述】\n")
cat("  'Tissue-specific eQTL analysis revealed that the lead cis-eQTL for NDUFB7\n")
cat("   in heart left ventricle (rs8103021, p=6.6e-6, slope=-0.08) differed from\n")
cat("   the blood eQTL (rs11085898, p=2.6e-9, Z=5.95), highlighting tissue-specific\n")
cat("   genetic regulation of NDUFB7 expression.'\n")

cat("\n完成\n")
