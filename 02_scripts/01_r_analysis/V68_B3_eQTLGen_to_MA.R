df <- read.csv("~/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/eqtlgen_NDUFB7_extracted.csv", stringsAsFactors=FALSE)
message("▶ eQTLGen行数: ", nrow(df))
# SMR .ma格式: SNP A1 A2 freq b se p n
# Zscore → b近似: b = Zscore, se = 1.0 (因为 Z = b/se, 所以若se=1, b=Z)
out <- data.frame(
  SNP = df$SNP,
  A1 = df$AssessedAllele,
  A2 = df$OtherAllele,
  freq = 0.5,
  b = df$Zscore,
  se = 1.0,
  p = df$Pvalue,
  n = df$NrSamples,
  stringsAsFactors = FALSE
)
write.table(out, "~/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/smr_input/eQTLGen_NDUFB7.ma", row.names=FALSE, quote=FALSE, sep=" ")
message("✅ eQTLGen .ma 已生成: ", nrow(out), " SNPs")
