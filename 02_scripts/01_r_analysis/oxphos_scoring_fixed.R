library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("OXPHOS评分整合（修复版）: GSE57338 (n=313)\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 1. OXPHOS基因集 ----------
oxphos_genes <- c(
  "NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8",
  "NDUFV1","NDUFV2","NDUFA9","NDUFA10","NDUFA13",
  "NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
  "SDHA","SDHB","SDHC","SDHD",
  "UQCRC1","UQCRC2","UQCRH","CYC1",
  "COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX7A1",
  "ATP5A1","ATP5B","ATP5C1","ATP5D","ATP5O",
  "PPARGC1A","SIRT3","TFAM","NRF1","NFE2L2"  # NRF2官方符号是NFE2L2
)

cat("OXPHOS基因集:", length(oxphos_genes), "个\n")

# ---------- 2. 读取平台文件（使用B家族成功的策略）----------
plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)

# 找基因注释列
anno_col <- names(plat)[ncol(plat)]  # 通常最后一列是mrna_assignment
cat("注释列:", anno_col, "\n")

# 部分匹配：在注释中搜索OXPHOS基因名（B家族成功策略）
hits <- plat[grepl(paste(paste0(" // ", oxphos_genes, " // "), collapse="|"), plat[[anno_col]], ignore.case=TRUE)]
cat("OXPHOS探针候选:", nrow(hits), "\n")

# 解析基因符号（从注释中提取 // GENE // 格式）
get_gene <- function(x) {
  m <- regmatches(x, gregexpr("// [A-Z0-9]+ //", x))
  if(length(m[[1]]) > 0) {
    gsub(" // ", "", m[[1]][1])
  } else NA
}

hits$GeneSymbol <- sapply(hits[[anno_col]], get_gene)
hits <- hits[!is.na(hits$GeneSymbol)]
hits <- hits[hits$GeneSymbol %in% oxphos_genes]

# 去重
hits <- hits[!duplicated(hits$GeneSymbol)]
cat("去重后OXPHOS基因:", nrow(hits), ":", paste(hits$GeneSymbol, collapse=", "), "\n")

# ---------- 3. 提取表达 ----------
gse57338 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")

probe_ids <- as.character(hits[[1]])
expr_ids <- as.character(gse57338[[1]])
matched <- expr_ids %in% probe_ids

cat("表达矩阵匹配:", sum(matched), "/", nrow(hits), "\n")

if(sum(matched) == 0) {
  cat("❌ 无匹配，尝试直接搜索探针ID...\n")
  # 备用：在表达矩阵行名中搜索基因名
  # 可能需要重新读取带header的平台文件
  quit(status=1)
}

oxphos_expr <- gse57338[matched, ]
oxphos_expr$GeneSymbol <- hits$GeneSymbol[match(as.character(oxphos_expr[[1]]), probe_ids)]

# ---------- 4. 计算OXPHOS评分 ----------
sample_cols <- setdiff(names(oxphos_expr), c(names(oxphos_expr)[1], "GeneSymbol"))
oxphos_mat <- as.matrix(oxphos_expr[, ..sample_cols])
rownames(oxphos_mat) <- oxphos_expr$GeneSymbol

# z-score标准化后取均值
oxphos_scores <- colMeans(scale(t(oxphos_mat)), na.rm=TRUE)
cat("\nOXPHOS评分范围:", round(min(oxphos_scores), 2), "~", round(max(oxphos_scores), 2), "\n")

# ---------- 5. 与NDUFB7相关性 ----------
ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
ndufb7_vals <- as.numeric(unlist(ndufb7[, -1]))

cor_result <- cor(oxphos_scores, ndufb7_vals, use="complete.obs")
cat("\n【NDUFB7 vs OXPHOS评分】\n")
cat("  Pearson r:", round(cor_result, 3), "\n")
cat("  解释:", ifelse(cor_result > 0.3, "强正相关", ifelse(cor_result > 0.1, "中等正相关", "弱相关")), "\n")

# ---------- 6. 病因分层 ----------
pheno <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_phenotype.txt")
ds <- pheno$disease.status

cat("\n【病因分层OXPHOS评分】\n")
for(g in unique(ds)) {
  s <- oxphos_scores[ds == g]
  cat(g, "(n=", length(s), "): 均值=", round(mean(s), 3), ", SD=", round(sd(s), 3), "\n")
}

# 统计
kw <- kruskal.test(oxphos_scores ~ ds)
cat("\nKruskal-Wallis: p=", format(kw$p.value, digits=3), "\n")

# ---------- 7. 保存 ----------
result <- data.table(Sample=names(oxphos_scores), OXPHOS_Score=oxphos_scores, NDUFB7=ndufb7_vals, Disease=ds)
fwrite(result, "~/Projects/NDUFB7_HF_2026_04_20/03_results/oxphos_scores_fixed.txt", sep="\t")
cat("\n✅ 已保存: oxphos_scores_fixed.txt\n")

cat("\n完成\n")
