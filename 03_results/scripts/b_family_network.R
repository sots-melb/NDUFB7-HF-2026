library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("NDUFB7 B亚家族网络构建\n")
cat("=" ,rep("=", 69), "\n", sep="")

# B亚家族成员（Complex I accessory subunits）
b_family <- c("NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", 
              "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11")

cat("B亚家族成员:", length(b_family), "\n")
cat("  ", paste(b_family, collapse=", "), "\n")

# 从GSE57338提取B家族全部成员的表达
gse57338_expr <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")
gse57338_pheno <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_phenotype.txt")

# 需要平台文件映射探针ID到基因名
# 简化：使用GPL11532平台文件
plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)

# 查找B家族探针
gene_col <- NULL
for(col in names(plat)) {
  if(any(grepl("gene|symbol|Gene|Symbol", col, ignore.case=TRUE))) {
    gene_col <- col
    break
  }
}

if(!is.null(gene_col)) {
  b_probes <- plat[grepl(paste(paste0("^", b_family, "$"), collapse="|"), plat[[gene_col]], ignore.case=TRUE)]
  cat("\nB家族探针匹配:", nrow(b_probes), "\n")
  print(b_probes[, c(names(plat)[1], gene_col), with=FALSE])
  
  # 提取表达
  probe_ids <- as.character(b_probes[[1]])
  expr_ids <- as.character(gse57338_expr[[1]])
  
  found <- expr_ids %in% probe_ids
  cat("在表达矩阵中匹配:", sum(found), "\n")
  
  if(sum(found) > 0) {
    b_expr <- gse57338_expr[found, ]
    fwrite(b_expr, "~/Projects/NDUFB7_HF_2026_04_20/03_results/b_family_expression.txt", sep="\t")
    cat("✅ B家族表达已保存\n")
    
    # 计算B家族内部相关性
    b_mat <- as.matrix(b_expr[, -1])
    rownames(b_mat) <- b_expr[[1]]
    
    # 简化：用探针ID代替基因名（后续映射）
    cor_mat <- cor(t(b_mat), use="pairwise.complete.obs")
    cat("\nB家族内部相关性矩阵:\n")
    print(round(cor_mat, 2))
    
    # 保存相关性
    fwrite(as.data.table(cor_mat, keep.rownames="Probe"), 
           "~/Projects/NDUFB7_HF_2026_04_20/03_results/b_family_correlation.txt", sep="\t")
  }
} else {
  cat("⚠️ 未找到基因名列\n")
}

cat("\n完成\n")
