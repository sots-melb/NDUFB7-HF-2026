library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("NDUFB7验证数据集批量提取 (GDS4772 + GSE116250 + GSE55296)\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ========== 1. GSE116250 RPKM提取 ==========
cat("\n【1/3】GSE116250 RPKM (RNA-seq, n=64)\n")
gse116250_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/GSE116250_rpkm.txt.gz"

if(file.exists(gse116250_file)) {
  # 先查看第一行（表头）
  header <- fread(cmd=paste("zcat", gse116250_file, "| head -1"), header=FALSE)
  cat("表头样本数:", ncol(header)-1, "\n")
  
  # 提取NDUFB7行
  ndufb7_line <- system(paste("zcat", gse116250_file, "| grep -i '^NDUFB7\\|^ENSG00000099795'"), intern=TRUE)
  
  if(length(ndufb7_line) > 0) {
    cat("✅ 找到NDUFB7 RPKM数据\n")
    # 解析行
    parts <- strsplit(ndufb7_line[1], "\t")[[1]]
    cat("基因ID:", parts[1], "\n")
    cat("样本数:", length(parts)-1, "\n")
    
    # 保存
    writeLines(ndufb7_line, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_rpkm_line.txt")
    
    # 读取完整矩阵获取样本名
    expr_116250 <- fread(cmd=paste("zcat", gse116250_file), header=TRUE)
    ndufb7_expr_116250 <- expr_116250[grepl("^NDUFB7|^ENSG00000099795", expr_116250[[1]], ignore.case=TRUE)]
    
    if(nrow(ndufb7_expr_116250) > 0) {
      vals <- as.numeric(unlist(ndufb7_expr_116250[, -1]))
      cat("NDUFB7 RPKM统计:\n")
      cat("  均值:", round(mean(vals, na.rm=TRUE), 2), "\n")
      cat("  中位数:", round(median(vals, na.rm=TRUE), 2), "\n")
      cat("  范围:", round(min(vals, na.rm=TRUE), 2), "-", round(max(vals, na.rm=TRUE), 2), "\n")
      fwrite(ndufb7_expr_116250, file="~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250/NDUFB7_rpkm_extracted.txt", sep="\t")
    }
  } else {
    cat("⚠️ 未找到NDUFB7，查看基因名格式...\n")
    system(paste("zcat", gse116250_file, "| head -5"))
  }
} else {
  cat("❌ GSE116250文件不存在\n")
}

# ========== 2. GDS4772 curated提取 ==========
cat("\n【2/3】GDS4772 curated (Affymetrix, n=29)\n")
gds_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/GDS4772.soft.gz"

if(file.exists(gds_file)) {
  # 解压并提取表达矩阵部分
  lines <- readLines(gzfile(gds_file))
  
  # 找dataset_table
  start_idx <- grep("^!dataset_table_begin", lines)
  end_idx <- grep("^!dataset_table_end", lines)
  
  if(length(start_idx) > 0 && length(end_idx) > 0) {
    cat("找到curated矩阵: 行", start_idx, "-", end_idx, "\n")
    matrix_lines <- lines[(start_idx + 1):(end_idx - 1)]
    
    # 写入临时文件
    tmp <- tempfile()
    writeLines(matrix_lines, tmp)
    gds_expr <- fread(tmp, header=TRUE)
    cat("矩阵维度:", nrow(gds_expr), "genes x", ncol(gds_expr)-1, "samples\n")
    
    # 查找NDUFB7 (GDS通常用Gene Symbol)
    ndufb7_gds <- gds_expr[grepl("^NDUFB7$", gds_expr[[1]], ignore.case=TRUE)]
    cat("NDUFB7匹配:", nrow(ndufb7_gds), "行\n")
    
    if(nrow(ndufb7_gds) > 0) {
      vals <- as.numeric(unlist(ndufb7_gds[, -1]))
      cat("NDUFB7表达统计:\n")
      cat("  均值:", round(mean(vals, na.rm=TRUE), 4), "\n")
      cat("  中位数:", round(median(vals, na.rm=TRUE), 4), "\n")
      fwrite(ndufb7_gds, file="~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/NDUFB7_GDS4772.txt", sep="\t")
    }
    
    # 提取表型 (subset信息)
    subset_lines <- lines[grep("^!subset_description", lines)]
    cat("\n表型分组:\n")
    print(subset_lines)
  } else {
    cat("⚠️ 未找到dataset_table标记\n")
  }
} else {
  cat("❌ GDS4772文件不存在\n")
}

# ========== 3. GSE55296 count提取 ==========
cat("\n【3/3】GSE55296 count_data (RNA-seq, n=36)\n")
gse55296_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/GSE55296_count_data.txt.gz"

if(file.exists(gse55296_file)) {
  header <- fread(cmd=paste("zcat", gse55296_file, "| head -1"), header=FALSE)
  cat("表头样本数:", ncol(header)-1, "\n")
  
  ndufb7_line <- system(paste("zcat", gse55296_file, "| grep -i '^NDUFB7\\|^ENSG00000099795'"), intern=TRUE)
  
  if(length(ndufb7_line) > 0) {
    cat("✅ 找到NDUFB7 count数据\n")
    parts <- strsplit(ndufb7_line[1], "\t")[[1]]
    cat("基因ID:", parts[1], "\n")
    
    writeLines(ndufb7_line, "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_count_line.txt")
    
    expr_55296 <- fread(cmd=paste("zcat", gse55296_file), header=TRUE)
    ndufb7_55296 <- expr_55296[grepl("^NDUFB7|^ENSG00000099795", expr_55296[[1]], ignore.case=TRUE)]
    
    if(nrow(ndufb7_55296) > 0) {
      vals <- as.numeric(unlist(ndufb7_55296[, -1]))
      cat("NDUFB7 count统计:\n")
      cat("  均值:", round(mean(vals, na.rm=TRUE), 2), "\n")
      cat("  中位数:", round(median(vals, na.rm=TRUE), 2), "\n")
      fwrite(ndufb7_55296, file="~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE55296/NDUFB7_count_extracted.txt", sep="\t")
    }
  } else {
    cat("⚠️ 未找到NDUFB7\n")
  }
} else {
  cat("❌ GSE55296文件不存在\n")
}

cat("\n完成\n")
