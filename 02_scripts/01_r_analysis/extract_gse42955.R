library(data.table)

out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE42955"

cat("=" ,rep("=", 59), "\n", sep="")
cat("GSE42955: 提取表达矩阵 + NDUFB7 + 表型\n")
cat("=" ,rep("=", 59), "\n", sep="")

# 1. 读取series_matrix
sm_file <- file.path(out_dir, "GSE42955_series_matrix.txt.gz")
lines <- readLines(gzfile(sm_file))
cat("总行数:", length(lines), "\n")

# 找到矩阵
start_idx <- grep("!series_matrix_table_begin", lines)
end_idx <- grep("!series_matrix_table_end", lines)
cat("矩阵行:", start_idx, "-", end_idx, "\n")

matrix_lines <- lines[(start_idx + 1):(end_idx - 1)]
tmp_file <- tempfile()
writeLines(matrix_lines, tmp_file)

expr <- fread(tmp_file, header=TRUE)
cat("表达矩阵:", nrow(expr), "genes x", ncol(expr)-1, "samples\n")
cat("ID列:", names(expr)[1], "\n")
cat("ID示例:", head(expr[[1]], 5), "\n")

# 保存完整矩阵
fwrite(expr, file=file.path(out_dir, "GSE42955_series_matrix_expr.txt"), sep="\t")
cat("✅ 表达矩阵已保存\n")

# 2. 提取表型
pheno_lines <- lines[grep("^!Sample_characteristics", lines)]
cat("\n表型行数:", length(pheno_lines), "\n")

pheno_list <- list()
for(pl in pheno_lines) {
  parts <- strsplit(pl, "\t")[[1]]
  if(length(parts) > 1) {
    first_val <- gsub('"', '', parts[2])
    if(grepl(":", first_val)) {
      key <- trimws(strsplit(first_val, ":")[[1]][1])
      vals <- sapply(parts[-1], function(x) {
        x <- gsub('"', '', x)
        if(grepl(":", x)) trimws(strsplit(x, ":")[[1]][2]) else x
      })
      pheno_list[[key]] <- vals
    }
  }
}

cat("表型字段:", names(pheno_list), "\n")
pheno_df <- as.data.frame(pheno_list)
fwrite(pheno_df, file=file.path(out_dir, "GSE42955_phenotype.txt"), sep="\t")
cat("✅ 表型已保存\n")

# 显示关键分组
if("disease state" %in% names(pheno_list)) {
  cat("\n疾病状态:\n"); print(table(pheno_list[["disease state"]]))
}
if("tissue" %in% names(pheno_list)) {
  cat("\n组织:\n"); print(table(pheno_list[["tissue"]]))
}

# 3. 查找NDUFB7（需要平台文件映射）
cat("\n【NDUFB7查找】\n")
cat("GSE42955使用GPL6244 (Affymetrix Human Gene 1.0 ST)\n")
cat("探针ID格式可能与GSE57338(GPL11532)不同\n")
cat("需要先下载GPL6244平台文件进行映射\n")
cat("临时方案: 直接搜索表达矩阵中的ID是否包含NDUFB7相关信息\n")

# 检查ID格式
ids <- as.character(expr[[1]])
cat("ID格式示例:", head(ids, 10), "\n")

# 如果是数字ID（如7892501），需要平台文件
# 如果是其他格式，尝试直接匹配
ndufb7_rows <- expr[grepl("NDUFB7|ndufb7|ENSG00000099795", ids, ignore.case=TRUE)]
cat("直接匹配NDUFB7:", nrow(ndufb7_rows), "行\n")

if(nrow(ndufb7_rows) > 0) {
  print(ndufb7_rows[[1]])
  fwrite(ndufb7_rows, file=file.path(out_dir, "NDUFB7_expression_GSE42955_direct.txt"), sep="\t")
} else {
  cat("⚠️ 需要GPL6244平台文件映射\n")
  cat("下载: https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/suppl/GPL6244-17930.txt\n")
}

cat("\n完成\n")
