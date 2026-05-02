library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("GSE290577 Xenium分析: NDUFB7单细胞空间验证\n")
cat("=" ,rep("=", 69), "\n", sep="")

rds_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/04_spatial_geo/GSE290577/GSE290577_heart_spatial.rds"

if(!file.exists(rds_file)) {
  cat("❌ RDS文件不存在\n")
  quit(status=1)
}

cat("读取RDS:", rds_file, "\n")
xenium <- readRDS(rds_file)

cat("RDS对象类型:", class(xenium), "\n")
cat("对象结构:\n")
print(str(xenium))

# 如果是Seurat对象
if(require("Seurat", quietly=TRUE) && inherits(xenium, "Seurat")) {
  cat("\n【Seurat对象分析】\n")
  cat("细胞数:", ncol(xenium), "\n")
  cat("基因数:", nrow(xenium), "\n")
  cat("Assays:", names(xenium@assays), "\n")
  
  # 查找NDUFB7
  ndufb7_pattern <- "NDUFB7|ENSG00000099795"
  ndufb7_found <- grep(ndufb7_pattern, rownames(xenium), value=TRUE, ignore.case=TRUE)
  cat("NDUFB7匹配:", length(ndufb7_found), ":", paste(ndufb7_found, collapse=", "), "\n")
  
  if(length(ndufb7_found) > 0) {
    # 提取NDUFB7表达
    ndufb7_expr <- FetchData(xenium, vars=ndufb7_found[1])
    cat("NDUFB7表达范围:", round(min(ndufb7_expr), 3), "~", round(max(ndufb7_expr), 3), "\n")
    
    # 按细胞类型分组（如果有注释）
    if("cell_type" %in% colnames(xenium@meta.data)) {
      ct <- xenium@meta.data$cell_type
      cat("\n【细胞类型NDUFB7表达】\n")
      for(c in unique(ct)) {
        vals <- ndufb7_expr[ct == c, 1]
        cat(c, "(n=", length(vals), "): 均值=", round(mean(vals), 3), "\n")
      }
    }
    
    # 保存
    fwrite(data.table(Cell=rownames(ndufb7_expr), NDUFB7=ndufb7_expr[,1]),
           "~/Projects/NDUFB7_HF_2026_04_20/03_results/xenium_ndufb7_expression.txt", sep="\t")
    cat("\n✅ Xenium NDUFB7表达已保存\n")
  }
  
} else if(inherits(xenium, "data.frame") || inherits(xenium, "matrix")) {
  cat("\n【数据框/矩阵分析】\n")
  cat("维度:", nrow(xenium), "x", ncol(xenium), "\n")
  cat("列名:", head(colnames(xenium), 10), "\n")
  cat("行名:", head(rownames(xenium), 10), "\n")
  
  # 查找NDUFB7
  ndufb7_found <- grep("NDUFB7|ENSG00000099795", rownames(xenium), value=TRUE, ignore.case=TRUE)
  cat("NDUFB7匹配:", paste(ndufb7_found, collapse=", "), "\n")
}

cat("\n完成\n")
