library(GEOquery)

cat("尝试获取GSE116250表达矩阵...\n")
gse <- getGEO("GSE116250", GSEMatrix = TRUE)

if(length(gse) > 0) {
  cat("获取到", length(gse), "个表达矩阵\n")
  expr <- exprs(gse[[1]])
  cat("表达矩阵维度:", nrow(expr), "genes x", ncol(expr), "samples\n")
  
  # 保存
  out_dir <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE116250"
  fwrite(as.data.table(expr, keep.rownames="ID"), 
         file=file.path(out_dir, "GSE116250_expression.txt"), sep="\t")
  cat("✅ 表达矩阵已保存\n")
  
  # 查找NDUFB7
  ndufb7_rows <- rownames(expr)[grepl("NDUFB7|ENSG00000099795", rownames(expr), ignore.case=TRUE)]
  cat("NDUFB7匹配:", length(ndufb7_rows), "\n")
  if(length(ndufb7_rows) > 0) {
    cat("ID:", ndufb7_rows, "\n")
    cat("表达值:", expr[ndufb7_rows, 1:5], "\n")
  }
} else {
  cat("❌ GEO未提供处理后的表达矩阵，需要从SRA下载原始数据\n")
}
