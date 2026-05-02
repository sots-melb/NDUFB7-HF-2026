library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("修复GSE116250 + GSE55296: 多格式基因名搜索\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 1. GSE116250: 尝试Ensembl ID ----------
cat("\n【GSE116250】尝试Ensembl ID搜索\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/GSE116250_rpkm.txt.gz"

all_data <- fread(cmd=paste("zcat", f), header=TRUE)
cat("总矩阵:", nrow(all_data), "x", ncol(all_data)-1, "\n")
cat("第一列名:", names(all_data)[1], "\n")
cat("前5个ID:", head(all_data[[1]], 5), "\n")

# 方法A: 搜索ENSG00000099795
hits <- all_data[grepl("ENSG00000099795", all_data[[1]], ignore.case=TRUE)]
cat("Ensembl匹配:", nrow(hits), "\n")

# 方法B: 搜索任何包含NDUFB7的行（任意列）
if(nrow(hits) == 0) {
  for(col in names(all_data)[1:min(3, ncol(all_data))]) {
    hits <- all_data[grepl("NDUFB7", all_data[[col]], ignore.case=TRUE)]
    if(nrow(hits) > 0) {
      cat("在列'", col, "'找到NDUFB7:", nrow(hits), "\n")
      break
    }
  }
}

if(nrow(hits) > 0) {
  vals <- as.numeric(unlist(hits[, -1]))
  vals <- vals[!is.na(vals)]
  cat("有效样本:", length(vals), "\n")
  cat("NDUFB7 RPKM: 均值=", round(mean(vals),2), ", 中位数=", round(median(vals),2), "\n")
  cat("范围:", round(min(vals),2), "-", round(max(vals),2), "\n")
  fwrite(hits, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_rpkm_fixed.txt", sep="\t")
} else {
  cat("⚠️ GSE116250中未找到NDUFB7\n")
  cat("  可能原因: 基因名使用不同ID体系\n")
  cat("  建议: 查看series_matrix样本表型，确认表达矩阵行名格式\n")
}

# ---------- 2. GSE55296: 尝试多种格式 ----------
cat("\n【GSE55296】尝试多种格式搜索\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/GSE55296_count_data.txt.gz"

all_data <- fread(cmd=paste("zcat", f), header=TRUE)
cat("总矩阵:", nrow(all_data), "x", ncol(all_data)-1, "\n")
cat("第一列名:", names(all_data)[1], "\n")
cat("前5个ID:", head(all_data[[1]], 5), "\n")

# 搜索NDUFB7
hits <- all_data[grepl("NDUFB7|ENSG00000099795", all_data[[1]], ignore.case=TRUE)]
cat("第一列匹配:", nrow(hits), "\n")

if(nrow(hits) == 0) {
  for(col in names(all_data)[1:min(3, ncol(all_data))]) {
    hits <- all_data[grepl("NDUFB7|ENSG00000099795", all_data[[col]], ignore.case=TRUE)]
    if(nrow(hits) > 0) {
      cat("在列'", col, "'找到:", nrow(hits), "\n")
      break
    }
  }
}

if(nrow(hits) > 0) {
  vals <- as.numeric(unlist(hits[, -1]))
  vals <- vals[!is.na(vals)]
  cat("有效样本:", length(vals), "\n")
  cat("NDUFB7 count: 均值=", round(mean(vals),2), ", 中位数=", round(median(vals),2), "\n")
  fwrite(hits, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_count_fixed.txt", sep="\t")
} else {
  cat("⚠️ GSE55296中未找到NDUFB7\n")
}

cat("\n完成\n")
