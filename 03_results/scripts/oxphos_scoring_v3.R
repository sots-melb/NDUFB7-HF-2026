library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("OXPHOS评分整合 V3: 使用B家族成功策略\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 1. 定义OXPHOS基因 ----------
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

# ---------- 2. 读取平台文件（B家族成功策略）----------
plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)
anno_col <- names(plat)[ncol(plat)]

# 策略：在注释列中直接搜索OXPHOS基因名（部分匹配）
# 格式示例: NM_004545 // NDUFB1 // NADH dehydrogenase...
cat("注释列:", anno_col, "\n")
cat("注释示例:", substr(plat[[anno_col]][10000], 1, 100), "...\n")

# 搜索所有OXPHOS基因（部分匹配，如 // NDUFB7 //）
search_pattern <- paste(paste0(" // ", oxphos_genes, " // "), collapse="|")
hits <- plat[grepl(search_pattern, plat[[anno_col]], ignore.case=TRUE)]
cat("OXPHOS探针匹配:", nrow(hits), "\n")

# 解析基因符号（从 // GENE // 提取）
get_symbol <- function(x) {
  # 匹配 // GENE // 模式
  m <- regmatches(x, gregexpr("// [A-Z][A-Za-z0-9]+ //", x))
  if(length(m[[1]]) > 0) {
    # 提取所有候选，选择在oxphos_genes中的
    candidates <- gsub(" // ", "", m[[1]])
    found <- candidates[candidates %in% oxphos_genes]
    if(length(found) > 0) return(found[1])
  }
  return(NA)
}

hits$GeneSymbol <- sapply(hits[[anno_col]], get_symbol)
hits <- hits[!is.na(hits$GeneSymbol)]
hits <- hits[!duplicated(hits$GeneSymbol)]
cat("去重后:", nrow(hits), ":", paste(hits$GeneSymbol, collapse=", "), "\n")

# ---------- 3. 提取表达 ----------
gse57338 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")
probe_ids <- as.character(hits[[1]])
expr_ids <- as.character(gse57338[[1]])
matched <- expr_ids %in% probe_ids
cat("表达矩阵匹配:", sum(matched), "/", nrow(hits), "\n")

if(sum(matched) == 0) {
  cat("❌ 仍无匹配，检查ID格式...\n")
  cat("探针ID示例:", head(probe_ids, 3), "\n")
  cat("表达ID示例:", head(expr_ids, 3), "\n")
  quit(status=1)
}

oxphos_expr <- gse57338[matched, ]
oxphos_expr$GeneSymbol <- hits$GeneSymbol[match(as.character(oxphos_expr[[1]]), probe_ids)]

# ---------- 4. 计算评分 ----------
sample_cols <- setdiff(names(oxphos_expr), c(names(oxphos_expr)[1], "GeneSymbol"))
oxphos_mat <- as.matrix(oxphos_expr[, ..sample_cols])
rownames(oxphos_mat) <- oxphos_expr$GeneSymbol

oxphos_scores <- colMeans(scale(t(oxphos_mat)), na.rm=TRUE)
cat("\nOXPHOS评分范围:", round(min(oxphos_scores), 2), "~", round(max(oxphos_scores), 2), "\n")

# ---------- 5. 与NDUFB7相关性 ----------
ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
ndufb7_vals <- as.numeric(unlist(ndufb7[, -1]))

cor_r <- cor(oxphos_scores, ndufb7_vals, use="complete.obs")
cat("\n【NDUFB7 vs OXPHOS】 r =", round(cor_r, 3), "\n")

# ---------- 6. 病因分层 ----------
pheno <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_phenotype.txt")
ds <- pheno$disease.status

cat("\n【病因分层】\n")
for(g in unique(ds)) {
  s <- oxphos_scores[ds == g]
  cat(g, "(n=", length(s), "): 均值=", round(mean(s), 3), "\n")
}

kw <- kruskal.test(oxphos_scores ~ ds)
cat("\nKW p =", format(kw$p.value, digits=3), "\n")

# 保存
result <- data.table(Sample=names(oxphos_scores), OXPHOS=oxphos_scores, NDUFB7=ndufb7_vals, Disease=ds)
fwrite(result, "~/Projects/NDUFB7_HF_2026_04_20/03_results/oxphos_scores_v3.txt", sep="\t")
cat("\n✅ 已保存\n")

cat("\n完成\n")
