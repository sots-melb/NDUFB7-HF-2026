[cd ~/Projects/NDUFB7_HF_2026_04_20/03_results/mr_results

# 检查eQTLGen文件列名和Zscore
head -1 NDUFB7_exposure_eQTLGen_strongIV.csv
awk -F',' 'NR==2 {for(i=1;i<=NF;i++) print i, $i}' NDUFB7_exposure_eQTLGen_strongIV.csv | head -5

# 如果Zscore为空，从beta/se反推
cat > fix_zscore.R << 'REOF'
library(data.table)
f <- "NDUFB7_exposure_eQTLGen_strongIV.csv"
d <- fread(f)
cat("列名:", names(d), "\n")
cat("Zscore前3行:", head(d$Zscore, 3), "\n")
if(all(is.na(d$Zscore) | d$Zscore == "")) {
  cat("Zscore为空，从beta/se反推...\n")
  d$Zscore <- d$beta / d$se
  fwrite(d, f)
  cat("✅ Zscore已修复并保存\n")
}
