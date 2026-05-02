library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("GSE116250: 病因分层 + NDUFB7差异分析\n")
cat("=" ,rep("=", 69), "\n", sep="")

# 读取series_matrix解析表型
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/GSE116250_series_matrix.txt.gz"
lines <- readLines(gzfile(f))

# 提取disease行
disease_lines <- lines[grep("^!Sample_characteristics_ch1.*disease", lines)]
cat("Disease特征行:", length(disease_lines), "\n")

# 解析每个样本的disease状态
# 格式: !Sample_characteristics_ch1	"disease: non-failing"	"disease: DCM" ...
disease_parts <- strsplit(disease_lines[1], "\t")[[1]]
disease_vals <- gsub('"', '', disease_parts[-1])
cat("样本病因分布:\n")
print(table(disease_vals))

# 读取NDUFB7表达
ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_rpkm_fixed.txt")
expr_vals <- as.numeric(unlist(ndufb7[, -1]))
expr_vals <- expr_vals[!is.na(expr_vals)]

cat("\n样本数匹配:", length(expr_vals), "表达值 vs", length(disease_vals), "表型\n")

# 确保顺序一致（假设series_matrix和表达矩阵样本顺序相同）
if(length(expr_vals) == length(disease_vals)) {
  df <- data.table(
    Expression = expr_vals,
    Disease = disease_vals
  )
  
  cat("\n【GSE116250病因分层】\n")
  for(d in unique(disease_vals)) {
    subset <- df[Disease == d]
    cat(d, "(n=", nrow(subset), "): ")
    cat("均值=", round(mean(subset$Expression), 2), ", ")
    cat("中位数=", round(median(subset$Expression), 2), "\n")
  }
  
  # 统计检验
  nf_vals <- df[Disease == "non-failing"]$Expression
  dcm_vals <- df[Disease == "dilated cardiomyopathy"]$Expression
  icm_vals <- df[Disease == "ischemic cardiomyopathy"]$Expression
  
  if(length(nf_vals) > 0 && length(dcm_vals) > 0) {
    wt_dcm <- wilcox.test(dcm_vals, nf_vals)
    cat("\nDCM vs NF: p=", format(wt_dcm$p.value, digits=4), "\n")
  }
  if(length(nf_vals) > 0 && length(icm_vals) > 0) {
    wt_icm <- wilcox.test(icm_vals, nf_vals)
    cat("ICM vs NF: p=", format(wt_icm$p.value, digits=4), "\n")
  }
  if(length(dcm_vals) > 0 && length(icm_vals) > 0) {
    wt_dcm_icm <- wilcox.test(dcm_vals, icm_vals)
    cat("DCM vs ICM: p=", format(wt_dcm_icm$p.value, digits=4), "\n")
  }
  if(length(unique(disease_vals)) >= 2) {
    kw <- kruskal.test(Expression ~ Disease, data=df)
    cat("Kruskal-Wallis: p=", format(kw$p.value, digits=4), "\n")
  }
  
  # 保存
  fwrite(df, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_stratified.txt", sep="\t")
  cat("\n✅ 分层数据已保存\n")
} else {
  cat("⚠️ 样本数不匹配，需要手动匹配样本顺序\n")
}

cat("\n完成\n")
