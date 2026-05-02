library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("GSE55296: 基于readme的精确病因分层\n")
cat("=" ,rep("=", 69), "\n", sep="")

# 读取readme构建映射
readme <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/GSE55296_processed_data_readme.txt", skip=1)
names(readme) <- c("GSM", "Title", "Column")
cat("Readme样本数:", nrow(readme), "\n")

# 从Title提取分组
readme$Group <- gsub(" rep[0-9]+", "", readme$Title)
cat("分组分布:\n")
print(table(readme$Group))

# 读取NDUFB7表达
ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_count_fixed.txt")
expr_vals <- as.numeric(unlist(ndufb7[, -1]))
expr_vals <- expr_vals[!is.na(expr_vals)]

# 读取count_data完整列名
count_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/GSE55296_count_data.txt.gz"
header <- system(paste("zcat", count_file, "| head -1"), intern=TRUE)
count_cols <- strsplit(header, "\t")[[1]][-1]
cat("\nCount data列数:", length(count_cols), "\n")

# 匹配readme列名到count data列名
# readme的Column列如"G114", "G12"等
# count data列名可能需要trim

matched <- readme[Column %in% trimws(count_cols)]
cat("匹配到count data的样本:", nrow(matched), "\n")

if(nrow(matched) > 0) {
  # 提取匹配样本的表达值
  matched$Expression <- NA
  for(i in 1:nrow(matched)) {
    col_idx <- which(trimws(count_cols) == matched$Column[i])
    if(length(col_idx) > 0) {
      matched$Expression[i] <- expr_vals[col_idx]
    }
  }
  
  cat("\n【GSE55296病因分层】\n")
  for(g in unique(matched$Group)) {
    subset <- matched[Group == g]
    vals <- subset$Expression[!is.na(subset$Expression)]
    cat(g, "(n=", length(vals), "): ")
    cat("均值=", round(mean(vals), 2), ", ")
    cat("中位数=", round(median(vals), 2), "\n")
  }
  
  # 统计检验
  ctl <- matched[Group == "Control"]$Expression
  dcm <- matched[Group == "Dilated"]$Expression
  icm <- matched[Group == "Ischemic"]$Expression
  
  cat("\n【统计比较】\n")
  if(length(ctl) > 0 && length(dcm) > 0) {
    wt <- wilcox.test(dcm, ctl)
    cat("DCM vs Control: p=", format(wt$p.value, digits=4), "\n")
  }
  if(length(ctl) > 0 && length(icm) > 0) {
    wt <- wilcox.test(icm, ctl)
    cat("Ischemic vs Control: p=", format(wt$p.value, digits=4), "\n")
  }
  if(length(dcm) > 0 && length(icm) > 0) {
    wt <- wilcox.test(dcm, icm)
    cat("DCM vs Ischemic: p=", format(wt$p.value, digits=4), "\n")
  }
  
  # 保存
  fwrite(matched, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_stratified.txt", sep="\t")
  cat("\n✅ 分层数据已保存\n")
} else {
  cat("⚠️ 列名匹配失败\n")
  cat("Readme列名示例:", head(readme$Column, 5), "\n")
  cat("Count列名示例:", head(trimws(count_cols), 5), "\n")
}

cat("\n完成\n")
