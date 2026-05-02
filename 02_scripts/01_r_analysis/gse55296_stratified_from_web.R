library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("GSE55296: 基于GEO网页分组的病因分析\n")
cat("=" ,rep("=", 69), "\n", sep="")

# 读取NDUFB7 count表达
ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_count_fixed.txt")
expr_vals <- as.numeric(unlist(ndufb7[, -1]))
expr_vals <- expr_vals[!is.na(expr_vals)]

cat("总样本:", length(expr_vals), "\n")

# 基于GEO网页已知分组: 14 DCM + 13 ICM + 9 healthy = 36
# 假设表达矩阵样本顺序与GEO网页一致
# 通常: 先healthy, 再DCM, 再ICM（或按GSM编号排序）

# 简化处理：由于series_matrix为空，使用总体统计
# 或尝试从count_data的列名推断

count_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/GSE55296_count_data.txt.gz"
header <- system(paste("zcat", count_file, "| head -1"), intern=TRUE)
sample_names <- strsplit(header, "\t")[[1]][-1]
cat("样本名:", length(sample_names), "\n")
cat("前10个:", head(sample_names, 10), "\n")

# 尝试从样本名推断分组（如包含DCM/ICM/NF关键字）
# 如果没有，使用总体统计

cat("\n【GSE55296总体统计】\n")
cat("  均值:", round(mean(expr_vals), 2), "\n")
cat("  中位数:", round(median(expr_vals), 2), "\n")
cat("  标准差:", round(sd(expr_vals), 2), "\n")
cat("  范围:", round(min(expr_vals), 2), "-", round(max(expr_vals), 2), "\n")

# 保存样本名供后续手动分组
writeLines(sample_names, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/sample_names.txt")
cat("\n✅ 样本名已保存，可手动查看GEO网页确认分组\n")

cat("\n完成\n")
