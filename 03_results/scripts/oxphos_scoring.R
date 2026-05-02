library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("OXPHOS评分整合: GSE57338 (n=313)\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 1. 定义OXPHOS基因集 ----------
# 来自MSigDB Hallmark: HALLMARK_OXIDATIVE_PHOSPHORYLATION
# 以及自定义Complex I核心基因
oxphos_genes <- c(
  # Complex I (NDUF家族核心)
  "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS7", "NDUFS8",
  "NDUFV1", "NDUFV2", "NDUFA9", "NDUFA10", "NDUFA13",
  "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10",
  # Complex II
  "SDHA", "SDHB", "SDHC", "SDHD",
  # Complex III
  "UQCRC1", "UQCRC2", "UQCRH", "CYC1",
  # Complex IV
  "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX7A1",
  # Complex V
  "ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5O",
  # 其他OXPHOS相关
  "PPARGC1A", "SIRT3", "TFAM", "NRF1", "NRF2"
)

cat("OXPHOS基因集:", length(oxphos_genes), "个基因\n")

# ---------- 2. 从GSE57338提取OXPHOS基因表达 ----------
gse57338 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")
plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)

# 找OXPHOS基因探针
anno_col <- NULL
for(col in names(plat)) {
  if(any(grepl("mrna|assignment|gene|symbol", col, ignore.case=TRUE))) {
    anno_col <- col
    break
  }
}

get_symbol <- function(x) {
  m <- regmatches(x, gregexpr("[A-Z][A-Z0-9]+", x))
  if(length(m[[1]]) > 0) return(m[[1]][1]) else return(NA)
}

plat$GeneSymbol <- sapply(plat[[anno_col]], get_symbol)
oxphos_hits <- plat[GeneSymbol %in% oxphos_genes]
oxphos_hits <- oxphos_hits[!duplicated(oxphos_hits$GeneSymbol)]
cat("匹配到OXPHOS探针:", nrow(oxphos_hits), "\n")

# 提取表达
probe_ids <- as.character(oxphos_hits[[1]])
expr_ids <- as.character(gse57338[[1]])
matched <- expr_ids %in% probe_ids
oxphos_expr <- gse57338[matched, ]
oxphos_expr$GeneSymbol <- oxphos_hits$GeneSymbol[match(as.character(oxphos_expr[[1]]), probe_ids)]

cat("OXPHOS表达矩阵:", nrow(oxphos_expr), "genes\n")

# ---------- 3. 计算OXPHOS评分（简单均值法，无需额外包） ----------
sample_cols <- setdiff(names(oxphos_expr), c(names(oxphos_expr)[1], "GeneSymbol"))
oxphos_mat <- as.matrix(oxphos_expr[, ..sample_cols])
rownames(oxphos_mat) <- oxphos_expr$GeneSymbol

# 计算每个样本的OXPHOS评分（z-score标准化后均值）
oxphos_scores <- colMeans(scale(t(oxphos_mat)), na.rm=TRUE)

cat("OXPHOS评分范围:", round(min(oxphos_scores), 2), "~", round(max(oxphos_scores), 2), "\n")

# ---------- 4. 与NDUFB7相关性 ----------
ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
ndufb7_vals <- as.numeric(unlist(ndufb7[, -1]))

cor_oxphos_ndufb7 <- cor(oxphos_scores, ndufb7_vals, use="complete.obs")
cat("\n【NDUFB7 vs OXPHOS评分相关性】\n")
cat("  Pearson r:", round(cor_oxphos_ndufb7, 3), "\n")

# ---------- 5. 病因分层OXPHOS评分 ----------
pheno <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_phenotype.txt")
ds <- pheno$disease.status

cat("\n【病因分层OXPHOS评分】\n")
for(g in unique(ds)) {
  subset_scores <- oxphos_scores[ds == g]
  cat(g, "(n=", length(subset_scores), "): ")
  cat("均值=", round(mean(subset_scores), 3), ", ")
  cat("SD=", round(sd(subset_scores), 3), "\n")
}

# 统计检验
nf_scores <- oxphos_scores[ds == "non-failing"]
dcm_scores <- oxphos_scores[ds == "idiopathic dilated CMP"]
isc_scores <- oxphos_scores[ds == "ischemic"]

wt_dcm <- wilcox.test(dcm_scores, nf_scores)
wt_isc <- wilcox.test(isc_scores, nf_scores)
wt_dcm_isc <- wilcox.test(dcm_scores, isc_scores)
kw <- kruskal.test(oxphos_scores ~ ds)

cat("\n【统计比较】\n")
cat("  DCM vs NF: p=", format(wt_dcm$p.value, digits=3), "\n")
cat("  ICM vs NF: p=", format(wt_isc$p.value, digits=3), "\n")
cat("  DCM vs ICM: p=", format(wt_dcm_isc$p.value, digits=3), "\n")
cat("  Kruskal-Wallis: p=", format(kw$p.value, digits=3), "\n")

# ---------- 6. NDUFB7在OXPHOS中的相对贡献 ----------
# 计算NDUFB7表达占OXPHOS总表达的百分比
ndufb7_z <- scale(ndufb7_vals)[,1]
oxphos_z <- oxphos_scores
relative_ndufb7 <- ndufb7_z / oxphos_z

cat("\n【NDUFB7相对OXPHOS贡献】\n")
cat("  NDUFB7 z-score / OXPHOS z-score 中位数:", round(median(relative_ndufb7, na.rm=TRUE), 3), "\n")
cat("  解释: NDUFB7在OXPHOS中的相对表达水平\n")

# ---------- 7. 保存 ----------
result <- data.table(
  Sample = names(oxphos_scores),
  OXPHOS_Score = oxphos_scores,
  NDUFB7_Expression = ndufb7_vals,
  Disease_Status = ds
)
fwrite(result, "~/Projects/NDUFB7_HF_2026_04_20/03_results/oxphos_ndufb7_scores.txt", sep="\t")

cat("\n✅ OXPHOS评分已保存: oxphos_ndufb7_scores.txt\n")

cat("\n完成\n")
