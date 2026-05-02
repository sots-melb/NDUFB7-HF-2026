library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("批量提取: GDS4772 + GSE116250 + GSE55296\n")
cat("=" ,rep("=", 69), "\n", sep="")

# --- 1. GSE116250 RPKM ---
cat("\n【1/3】GSE116250 RPKM\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/GSE116250_rpkm.txt.gz"
lines <- system(paste("zcat", f, "| grep -i '^NDUFB7\\|^ENSG00000099795'"), intern=TRUE)
if(length(lines) > 0) {
  parts <- strsplit(lines[1], "\t")[[1]]
  vals <- as.numeric(parts[-1])
  cat("样本数:", length(vals), "\n")
  cat("NDUFB7 RPKM: 均值=", round(mean(vals),2), ", 中位数=", round(median(vals),2), "\n")
  
  # 保存简洁版
  writeLines(lines, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_rpkm_extracted.txt")
} else cat("未找到\n")

# --- 2. GDS4772 ---
cat("\n【2/3】GDS4772 curated\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/GDS4772.soft.gz"
lines <- readLines(gzfile(f))
start <- grep("^!dataset_table_begin", lines)
end <- grep("^!dataset_table_end", lines)
if(length(start)>0) {
  mat <- lines[(start+1):(end-1)]
  tmp <- tempfile(); writeLines(mat, tmp)
  gds <- fread(tmp, header=TRUE)
  cat("矩阵:", nrow(gds), "x", ncol(gds)-1, "\n")
  
  # 找NDUFB7
  idcol <- names(gds)[1]
  hits <- gds[grepl("^NDUFB7$", gds[[idcol]], ignore.case=TRUE)]
  if(nrow(hits)==0 && ncol(gds)>1) {
    hits <- gds[grepl("NDUFB7", gds[[2]], ignore.case=TRUE)]
  }
  cat("NDUFB7匹配:", nrow(hits), "\n")
  if(nrow(hits)>0) {
    vals <- as.numeric(unlist(hits[, -1]))
    cat("表达: 均值=", round(mean(vals),4), ", 中位数=", round(median(vals),4), "\n")
    fwrite(hits, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_GDS4772.txt", sep="\t")
  }
}

# --- 3. GSE55296 count ---
cat("\n【3/3】GSE55296 count\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/GSE55296_count_data.txt.gz"
lines <- system(paste("zcat", f, "| grep -i '^NDUFB7\\|^ENSG00000099795'"), intern=TRUE)
if(length(lines)>0) {
  parts <- strsplit(lines[1], "\t")[[1]]
  vals <- as.numeric(parts[-1])
  cat("样本数:", length(vals), "\n")
  cat("NDUFB7 count: 均值=", round(mean(vals),2), ", 中位数=", round(median(vals),2), "\n")
  writeLines(lines, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_count_extracted.txt")
} else cat("未找到\n")

cat("\n完成\n")
