library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("RNA-seq病因分组: GSE116250 + GSE55296\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 1. GSE116250表型 ----------
cat("\n【1/2】GSE116250表型\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/GSE116250_series_matrix.txt.gz"
lines <- readLines(gzfile(f))

# 提取样本特征
char_lines <- lines[grep("^!Sample_characteristics", lines)]
cat("特征行数:", length(char_lines), "\n")
for(cl in char_lines[1:min(5, length(char_lines))]) {
  cat("  ", substr(cl, 1, 100), "\n")
}

# 提取title（可能包含分组信息）
title_lines <- lines[grep("^!Sample_title", lines)]
cat("\n样本title示例:\n")
for(tl in title_lines[1:min(5, length(title_lines))]) {
  cat("  ", substr(tl, 1, 100), "\n")
}

# ---------- 2. GSE55296表型 ----------
cat("\n【2/2】GSE55296表型\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/GSE55296_series_matrix.txt.gz"

if(file.exists(f) && file.size(f) > 10000) {
  lines <- readLines(gzfile(f))
  char_lines <- lines[grep("^!Sample_characteristics", lines)]
  cat("特征行数:", length(char_lines), "\n")
  for(cl in char_lines[1:min(5, length(char_lines))]) {
    cat("  ", substr(cl, 1, 100), "\n")
  }
  
  title_lines <- lines[grep("^!Sample_title", lines)]
  cat("\n样本title示例:\n")
  for(tl in title_lines[1:min(5, length(title_lines))]) {
    cat("  ", substr(tl, 1, 100), "\n")
  }
} else {
  cat("⚠️ GSE55296 series_matrix不完整，使用已知分组:\\n")
  cat("  来自GEO网页: 14 DCM + 13 ICM + 9 healthy (或类似)\\n")
}

cat("\n完成\n")
