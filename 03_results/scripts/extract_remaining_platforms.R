library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("批量提取: GDS4772 + GSE116250 + GSE55296\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ========== 1. GSE116250 RPKM ==========
cat("\n【1/3】GSE116250 RNA-seq RPKM (n=64)\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/GSE116250_rpkm.txt.gz"

# 读取表头获取样本名
header <- system(paste("zcat", f, "| head -1"), intern=TRUE)
sample_names <- strsplit(header, "\t")[[1]][-1]
cat("样本数:", length(sample_names), "\n")

# 提取NDUFB7
cmd <- paste("zcat", f, "| grep -i '^NDUFB7\\|^ENSG00000099795'")
ndufb7_line <- system(cmd, intern=TRUE)

if(length(ndufb7_line) > 0) {
  parts <- strsplit(ndufb7_line[1], "\t")[[1]]
  gene_id <- parts[1]
  vals <- as.numeric(parts[-1])
  
  cat("✅ NDUFB7 found:", gene_id, "\n")
  cat("  均值:", round(mean(vals, na.rm=TRUE), 2), "\n")
  cat("  中位数:", round(median(vals, na.rm=TRUE), 2), "\n")
  cat("  范围:", round(min(vals, na.rm=TRUE), 2), "-", round(max(vals, na.rm=TRUE), 2), "\n")
  
  # 保存
  df <- data.table(ID=gene_id, t(vals))
  setnames(df, c("ID", sample_names))
  fwrite(df, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_rpkm_extracted.txt", sep="\t")
  
  # 保存简洁版供整合
  writeLines(paste(gene_id, paste(vals, collapse="\t"), sep="\t"), 
             "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_rpkm_line.txt")
} else {
  cat("⚠️ NDUFB7 not found in RPKM\n")
}

# ========== 2. GDS4772 curated ==========
cat("\n【2/3】GDS4772 Affymetrix curated (n=29)\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/GDS4772.soft.gz"

lines <- readLines(gzfile(f))
start_idx <- grep("^!dataset_table_begin", lines)
end_idx <- grep("^!dataset_table_end", lines)

if(length(start_idx) > 0) {
  mat_lines <- lines[(start_idx+1):(end_idx-1)]
  tmp <- tempfile(); writeLines(mat_lines, tmp)
  gds <- fread(tmp, header=TRUE)
  cat("矩阵:", nrow(gds), "x", ncol(gds)-1, "\n")
  
  # 查找NDUFB7 (GDS通常用Gene Symbol作为第一列或ID_REF)
  id_col <- names(gds)[1]
  ndufb7_rows <- gds[grepl("^NDUFB7$", gds[[id_col]], ignore.case=TRUE)]
  
  if(nrow(ndufb7_rows) == 0) {
    # 尝试第二列
    if(ncol(gds) > 1) {
      id_col2 <- names(gds)[2]
      ndufb7_rows <- gds[grepl("NDUFB7", gds[[id_col2]], ignore.case=TRUE)]
    }
  }
  
  cat("NDUFB7匹配:", nrow(ndufb7_rows), "\n")
  if(nrow(ndufb7_rows) > 0) {
    vals <- as.numeric(unlist(ndufb7_rows[, -1]))
    cat("  均值:", round(mean(vals, na.rm=TRUE), 4), "\n")
    cat("  中位数:", round(median(vals, na.rm=TRUE), 4), "\n")
    fwrite(ndufb7_rows, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_GDS4772.txt", sep="\t")
  }
  
  # 提取subset表型
  subset_desc <- lines[grep("^!subset_description", lines)]
  subset_type <- lines[grep("^!subset_type", lines)]
  cat("\n表型分组:\n")
  for(i in 1:length(subset_desc)) {
    cat(" ", subset_desc[i], "\n")
  }
} else {
  cat("⚠️ No dataset_table found\n")
}

# ========== 3. GSE55296 count_data ==========
cat("\n【3/3】GSE55296 RNA-seq count (n=36)\n")
f <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/GSE55296_count_data.txt.gz"

header <- system(paste("zcat", f, "| head -1"), intern=TRUE)
sample_names_55296 <- strsplit(header, "\t")[[1]][-1]
cat("样本数:", length(sample_names_55296), "\n")

cmd <- paste("zcat", f, "| grep -i '^NDUFB7\\|^ENSG00000099795'")
ndufb7_line <- system(cmd, intern=TRUE)

if(length(ndufb7_line) > 0) {
  parts <- strsplit(ndufb7_line[1], "\t")[[1]]
  vals <- as.numeric(parts[-1])
  cat("✅ NDUFB7 found:", parts[1], "\n")
  cat("  均值:", round(mean(vals, na.rm=TRUE), 2), "\n")
  cat("  中位数:", round(median(vals, na.rm=TRUE), 2), "\n")
  
  df <- data.table(ID=parts[1], t(vals))
  setnames(df, c("ID", sample_names_55296))
  fwrite(df, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_count_extracted.txt", sep="\t")
} else {
  cat("⚠️ NDUFB7 not found in counts\n")
}

cat("\n完成\n")
