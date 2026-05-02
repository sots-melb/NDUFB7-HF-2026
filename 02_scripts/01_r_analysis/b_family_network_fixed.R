library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("B亚家族网络构建（修复版）\n")
cat("=" ,rep("=", 69), "\n", sep="")

b_family <- c("NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", 
              "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11")

cat("B亚家族:", paste(b_family, collapse=", "), "\n")

# 读取平台文件
plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)

cat("平台列名:", names(plat)[1:5], "...\n")
cat("平台维度:", nrow(plat), "x", ncol(plat), "\n")

# 查看基因注释列的结构（前5行）
cat("\n基因注释列样本:\n")
for(col in names(plat)[1:min(5, ncol(plat))]) {
  cat("  ", col, ":", substr(as.character(plat[[col]][1]), 1, 80), "\n")
}

# 策略：在所有列中搜索NDUFB（宽松匹配）
found_genes <- list()
for(col in names(plat)) {
  vals <- as.character(plat[[col]])
  for(gene in b_family) {
    hits <- grep(gene, vals, ignore.case=TRUE)
    if(length(hits) > 0) {
      if(is.null(found_genes[[gene]])) {
        found_genes[[gene]] <- data.table(
          Gene = gene,
          Column = col,
          Row = hits[1],
          ProbeID = as.character(plat[[1]][hits[1]]),
          Annotation = substr(vals[hits[1]], 1, 100)
        )
      }
    }
  }
}

if(length(found_genes) > 0) {
  b_map <- rbindlist(found_genes)
  cat("\nB家族探针匹配:\n")
  print(b_map)
  
  # 去重（每个基因取第一个匹配）
  b_unique <- b_map[!duplicated(Gene)]
  cat("\n去重后:", nrow(b_unique), "基因\n")
  
  # 提取GSE57338表达
  gse57338 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")
  
  probe_ids <- b_unique$ProbeID
  expr_ids <- as.character(gse57338[[1]])
  
  found_mask <- expr_ids %in% probe_ids
  cat("表达矩阵匹配:", sum(found_mask), "/", nrow(b_unique), "\n")
  
  if(sum(found_mask) > 0) {
    b_expr <- gse57338[found_mask, ]
    b_expr$GeneSymbol <- b_unique$Gene[match(as.character(b_expr[[1]]), b_unique$ProbeID)]
    
    fwrite(b_expr, "~/Projects/NDUFB7_HF_2026_04_20/03_results/b_family_expression.txt", sep="\t")
    cat("✅ B家族表达已保存\n")
    
    # 计算相关性
    b_mat <- as.matrix(b_expr[, -c(1, ncol(b_expr))])
    rownames(b_mat) <- b_expr$GeneSymbol
    
    b_mat_complete <- b_mat[complete.cases(b_mat), ]
    cat("完整基因:", nrow(b_mat_complete), "\n")
    
    if(nrow(b_mat_complete) >= 3) {
      cor_mat <- cor(t(b_mat_complete), use="pairwise.complete.obs")
      cat("\n相关性矩阵:\n")
      print(round(cor_mat, 2))
      
      fwrite(as.data.table(cor_mat, keep.rownames="Gene"), 
             "~/Projects/NDUFB7_HF_2026_04_20/03_results/b_family_correlation.txt", sep="\t")
      
      if("NDUFB7" %in% rownames(cor_mat)) {
        ndufb7_corr <- sort(cor_mat["NDUFB7", ], decreasing=TRUE)
        cat("\n【NDUFB7共表达排名】\n")
        print(round(ndufb7_corr, 3))
      }
    }
  }
} else {
  cat("⚠️ 未找到B家族探针\n")
}

cat("\n完成\n")
