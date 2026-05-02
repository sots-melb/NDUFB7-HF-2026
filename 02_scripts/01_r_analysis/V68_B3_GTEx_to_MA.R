df <- read.csv("~/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/gtex_NDUFB7_HLV_v11.csv", stringsAsFactors=FALSE)
message("▶ GTEx行数: ", nrow(df))
# variant_id: chr19_14061916_C_T_b38 → 提取等位基因
# 假设第二个等位基因是效应等位基因（替代等位基因）
alleles <- do.call(rbind, strsplit(df$variant_id, "_"))
# chr, pos, A1, A2, build
out <- data.frame(
  SNP = paste0(alleles[,1], ":", alleles[,2]),  # chr:pos
  A1 = alleles[,4],  # 替代等位基因（效应）
  A2 = alleles[,3],  # 参考等位基因
  freq = df$af,
  b = df$slope,
  se = df$slope_se,
  p = df$pval_nominal,
  n = df$ma_count,
  stringsAsFactors = FALSE
)
write.table(out, "~/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/smr_input/GTEx_HeartLV_NDUFB7.ma", row.names=FALSE, quote=FALSE, sep=" ")
message("✅ GTEx .ma 已生成: ", nrow(out), " SNPs")
