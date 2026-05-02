library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("修复NA: 四平台NDUFB7重新提取\n")
cat("=" ,rep("=", 69), "\n", sep="")

# --- 1. GSE116250 RPKM修复 ---
cat("\n【1/3】GSE116250 RPKM修复\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/GSE116250_rpkm.txt.gz"

# 方法: 用fread直接读取，指定基因名
cmd <- paste("zcat", f)
all_data <- fread(cmd=cmd, header=TRUE)
cat("总矩阵:", nrow(all_data), "x", ncol(all_data)-1, "\n")

# 查找NDUFB7（第一列可能是GeneSymbol或Ensembl）
id_col <- names(all_data)[1]
cat("ID列:", id_col, "\n")

# 方法A: 精确匹配NDUFB7
hits <- all_data[grepl("^NDUFB7$", all_data[[id_col]], ignore.case=TRUE)]
cat("精确匹配NDUFB7:", nrow(hits), "\n")

if(nrow(hits) == 0) {
  # 方法B: 包含匹配
  hits <- all_data[grepl("NDUFB7", all_data[[id_col]], ignore.case=TRUE)]
  cat("包含匹配:", nrow(hits), "\n")
}

if(nrow(hits) > 0) {
  vals <- as.numeric(unlist(hits[, -1]))
  vals <- vals[!is.na(vals)]
  cat("有效样本数:", length(vals), "\n")
  cat("NDUFB7 RPKM: 均值=", round(mean(vals),2), ", 中位数=", round(median(vals),2), "\n")
  cat("范围:", round(min(vals),2), "-", round(max(vals),2), "\n")
  
  # 保存
  fwrite(hits, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_rpkm_fixed.txt", sep="\t")
}

# --- 2. GDS4772修复 ---
cat("\n【2/3】GDS4772修复\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/GDS4772.soft.gz"

lines <- readLines(gzfile(f))
start <- grep("^!dataset_table_begin", lines)
end <- grep("^!dataset_table_end", lines)

mat_lines <- lines[(start+1):(end-1)]
tmp <- tempfile()
writeLines(mat_lines, tmp)

gds <- fread(tmp, header=TRUE)
cat("矩阵:", nrow(gds), "x", ncol(gds)-1, "\n")
cat("列名:", names(gds)[1:5], "...\n")

# GDS4772第一列是ID_REF（数字探针ID），需要找基因名列
# 通常第二列或最后一列有基因注释
for(i in 1:min(3, ncol(gds))) {
  col <- names(gds)[i]
  hits <- gds[grepl("NDUFB7", gds[[col]], ignore.case=TRUE)]
  if(nrow(hits) > 0) {
    cat("在第", i, "列('", col, "')找到NDUFB7:", nrow(hits), "\n")
    # 提取数值（跳过ID列和注释列）
    val_cols <- setdiff(names(gds), c(col, "ID_REF"))[1:(ncol(gds)-2)]
    vals <- as.numeric(unlist(hits[, ..val_cols]))
    vals <- vals[!is.na(vals)]
    cat("  有效样本:", length(vals), "\n")
    cat("  均值:", round(mean(vals),4), ", 中位数:", round(median(vals),4), "\n")
    fwrite(hits, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_GDS4772_fixed.txt", sep="\t")
    break
  }
}

# --- 3. GSE55296 count修复 ---
cat("\n【3/3】GSE55296 count修复\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/GSE55296_count_data.txt.gz"

all_data <- fread(cmd=paste("zcat", f), header=TRUE)
cat("总矩阵:", nrow(all_data), "x", ncol(all_data)-1, "\n")

id_col <- names(all_data)[1]
hits <- all_data[grepl("^NDUFB7$", all_data[[id_col]], ignore.case=TRUE)]

if(nrow(hits) == 0) {
  hits <- all_data[grepl("NDUFB7", all_data[[id_col]], ignore.case=TRUE)]
}

cat("NDUFB7匹配:", nrow(hits), "\n")
if(nrow(hits) > 0) {
  vals <- as.numeric(unlist(hits[, -1]))
  vals <- vals[!is.na(vals)]
  cat("有效样本:", length(vals), "\n")
  cat("NDUFB7 count: 均值=", round(mean(vals),2), ", 中位数=", round(median(vals),2), "\n")
  fwrite(hits, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_count_fixed.txt", sep="\t")
}

cat("\n完成\n")
