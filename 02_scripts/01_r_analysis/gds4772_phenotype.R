library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("GDS4772表型提取 + 病因分层\n")
cat("=" ,rep("=", 69), "\n", sep="")

f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/GDS4772.soft.gz"
lines <- readLines(gzfile(f))

# 提取subset信息（病因分组）
subset_desc <- lines[grep("^!subset_description", lines)]
subset_type <- lines[grep("^!subset_type", lines)]
subset_sample <- lines[grep("^!subset_sample_id", lines)]

cat("【GDS4772分组信息】\n")
for(i in 1:length(subset_desc)) {
  cat("\n--- Subset", i, "---\n")
  cat("  Description:", subset_desc[i], "\n")
  cat("  Type:", subset_type[i], "\n")
  cat("  Sample IDs:", substr(subset_sample[i], 1, 100), "...\n")
}

# 提取NDUFB7并按分组分析
ndufb7_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_GDS4772_fixed.txt"
if(file.exists(ndufb7_file)) {
  ndufb7 <- fread(ndufb7_file)
  cat("\n【NDUFB7表达值】\n")
  cat("列名:", names(ndufb7), "\n")
  
  # 显示所有样本的值
  val_cols <- setdiff(names(ndufb7), c("ID_REF", " IDENTIFIER "))
  cat("样本列数:", length(val_cols), "\n")
  
  vals <- as.numeric(unlist(ndufb7[, ..val_cols]))
  cat("表达值:", round(vals, 4), "\n")
}

cat("\n完成\n")
