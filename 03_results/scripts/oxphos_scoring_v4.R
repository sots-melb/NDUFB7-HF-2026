library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("OXPHOS评分整合 V4: 全列搜索策略（B家族成功方法）\n")
cat("=" ,rep("=", 69), "\n", sep="")

oxphos_genes <- c(
  "NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8",
  "NDUFV1","NDUFV2","NDUFA9","NDUFA10","NDUFA13",
  "NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
  "SDHA","SDHB","SDHC","SDHD",
  "UQCRC1","UQCRC2","UQCRH","CYC1",
  "COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX7A1",
  "ATP5A1","ATP5B","ATP5C1","ATP5D","ATP5O",
  "PPARGC1A","SIRT3","TFAM","NRF1","NFE2L2"
)

cat("OXPHOS基因:", length(oxphos_genes), "\n")

# 读取平台文件
plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)
cat("平台维度:", nrow(plat), "x", ncol(plat), "\n")
cat("列名:", names(plat)[1:5], "...\n")

# 策略：在所有列中搜索OXPHOS基因（B家族成功方法）
found_genes <- list()
for(col in names(plat)) {
  vals <- as.character(plat[[col]])
  for(gene in oxphos_genes) {
    hits <- grep(paste0(" // ", gene, " // "), vals, ignore.case=TRUE)
    if(length(hits) > 0) {
      found_genes[[gene]] <- data.table(
        Gene = gene,
        Column = col,
        Row = hits[1],
        ProbeID = as.character(plat[[1]][hits[1]]),
        Annotation = substr(vals[hits[1]], 1, 100)
      )
      break  # 找到即停止，避免重复
    }
  }
}

if(length(found_genes) == 0) {
  cat("⚠️ 严格模式无匹配，尝试宽松模式...\n")
  for(col in names(plat)) {
    vals <- as.character(plat[[col]])
    for(gene in oxphos_genes) {
      hits <- grep(gene, vals, ignore.case=TRUE)
      if(length(hits) > 0) {
        found_genes[[gene]] <- data.table(
          Gene = gene,
          Column = col,
          Row = hits[1],
          ProbeID = as.character(plat[[1]][hits[1]]),
          Annotation = substr(vals[hits[1]], 1, 100)
        )
        break
      }
    }
  }
}

cat("找到OXPHOS基因:", length(found_genes), "\n")
if(length(found_genes) > 0) {
  oxphos_df <- rbindlist(found_genes)
  print(oxphos_df[, .(Gene, ProbeID, Annotation)])
  
  # 提取表达
  gse57338 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")
  probe_ids <- oxphos_df$ProbeID
  expr_ids <- as.character(gse57338[[1]])
  matched <- expr_ids %in% probe_ids
  cat("\n表达矩阵匹配:", sum(matched), "/", length(probe_ids), "\n")
  
  if(sum(matched) > 0) {
    oxphos_expr <- gse57338[matched, ]
    sample_cols <- setdiff(names(oxphos_expr), names(oxphos_expr)[1])
    oxphos_mat <- as.matrix(oxphos_expr[, ..sample_cols])
    
    oxphos_scores <- colMeans(scale(t(oxphos_mat)), na.rm=TRUE)
    cat("OXPHOS评分范围:", round(min(oxphos_scores), 2), "~", round(max(oxphos_scores), 2), "\n")
    
    # NDUFB7相关性
    ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
    ndufb7_vals <- as.numeric(unlist(ndufb7[, -1]))
    cor_r <- cor(oxphos_scores, ndufb7_vals, use="complete.obs")
    cat("NDUFB7 vs OXPHOS r =", round(cor_r, 3), "\n")
    
    # 保存
    result <- data.table(Sample=names(oxphos_scores), OXPHOS=oxphos_scores, NDUFB7=ndufb7_vals)
    fwrite(result, "~/Projects/NDUFB7_HF_2026_04_20/03_results/oxphos_scores_v4.txt", sep="\t")
    cat("✅ 已保存\n")
  }
} else {
  cat("❌ 未找到任何OXPHOS基因\n")
}

cat("\n完成\n")
