library(data.table)
library(ggplot2)
library(reshape2)

cat("=" ,rep("=", 69), "\n", sep="")
cat("B家族 + 参照基因网络构建（Figure 3核心）\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 1. 定义基因列表 ----------
b_family <- c("NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", 
              "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11")

# 参照基因：核心亚基 + A家族 + 黄素蛋白 + 组装因子
ref_genes <- c("NDUFS4", "NDUFA9", "NDUFV1", "NDUFAF2", "NDUFC2")

all_genes <- c(b_family, ref_genes)
cat("目标基因:", length(all_genes), "个\n")
cat("  B家族:", paste(b_family, collapse=", "), "\n")
cat("  参照:", paste(ref_genes, collapse=", "), "\n")

# ---------- 2. 从GPL11532提取探针 ----------
plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)

# 找基因注释列（通常最后一列或包含mrna_assignment）
anno_col <- NULL
for(col in names(plat)) {
  if(any(grepl("mrna|assignment|gene|symbol", col, ignore.case=TRUE))) {
    anno_col <- col
    break
  }
}
if(is.null(anno_col)) anno_col <- names(plat)[ncol(plat)]
cat("注释列:", anno_col, "\n")

# 提取所有目标基因（部分匹配）
hits <- plat[grepl(paste(paste0("NDUF[A-Z]?[0-9]+"), collapse="|"), plat[[anno_col]], ignore.case=TRUE)]
cat("NDUF家族探针候选:", nrow(hits), "\n")

# 解析基因符号
get_symbol <- function(x) {
  m <- regmatches(x, gregexpr("NDUF[A-Z]?[0-9]+", x, ignore.case=TRUE))
  if(length(m[[1]]) > 0) return(toupper(m[[1]][1])) else return(NA)
}
hits$GeneSymbol <- sapply(hits[[anno_col]], get_symbol)
hits <- hits[!is.na(hits$GeneSymbol)]
hits <- hits[hits$GeneSymbol %in% all_genes]

# 去重：每个基因保留第一个探针
hits <- hits[!duplicated(hits$GeneSymbol)]
cat("去重后匹配:", nrow(hits), "/", length(all_genes), "\n")
print(hits[, c(names(plat)[1], anno_col, "GeneSymbol"), with=FALSE])

# ---------- 3. 提取GSE57338表达 ----------
gse57338 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")

probe_ids <- as.character(hits[[1]])
expr_ids <- as.character(gse57338[[1]])

matched <- expr_ids %in% probe_ids
cat("表达矩阵匹配:", sum(matched), "/", length(probe_ids), "\n")

if(sum(matched) == 0) {
  cat("❌ 无匹配，检查ID格式...\n")
  cat("探针ID示例:", head(probe_ids, 3), "\n")
  cat("表达ID示例:", head(expr_ids, 3), "\n")
  quit(status=1)
}

expr <- gse57338[matched, ]
expr$GeneSymbol <- hits$GeneSymbol[match(as.character(expr[[1]]), probe_ids)]

cat("成功提取基因:", paste(expr$GeneSymbol, collapse=", "), "\n")

# ---------- 4. 构建表达矩阵 ----------
# 明确指定样本列（第2列到倒数第2列，如果GeneSymbol是最后一列）
sample_col_names <- setdiff(names(expr), c(names(expr)[1], "GeneSymbol"))
cat("样本列数:", length(sample_col_names), "\n")

expr_mat <- as.matrix(expr[, ..sample_col_names])
rownames(expr_mat) <- expr$GeneSymbol
cat("表达矩阵维度:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

# 检查NA
na_per_gene <- rowSums(is.na(expr_mat))
cat("每基因NA数:\n")
print(na_per_gene)

# 只保留无NA的基因
keep_genes <- names(na_per_gene)[na_per_gene == 0]
expr_mat_clean <- expr_mat[keep_genes, , drop=FALSE]
cat("完整基因:", nrow(expr_mat_clean), ":", paste(rownames(expr_mat_clean), collapse=", "), "\n")

# ---------- 5. 计算相关性 ----------
cor_mat <- cor(t(expr_mat_clean), use="pairwise.complete.obs")
cat("\n相关性矩阵:\n")
print(round(cor_mat, 2))

# 保存
fwrite(as.data.table(cor_mat, keep.rownames="Gene"), 
       "~/Projects/NDUFB7_HF_2026_04_20/03_results/nduf_family_correlation_matrix.txt", sep="\t")
cat("\n✅ 相关性矩阵已保存\n")

# ---------- 6. NDUFB7特异性分析 ----------
if("NDUFB7" %in% rownames(cor_mat)) {
  cat("\n【NDUFB7共表达谱】\n")
  ndufb7_r <- sort(cor_mat["NDUFB7", ], decreasing=TRUE)
  print(round(ndufb7_r, 3))
  
  b_only <- ndufb7_r[grepl("^NDUFB[0-9]+$", names(ndufb7_r))]
  ref_only <- ndufb7_r[!grepl("^NDUFB[0-9]+$", names(ndufb7_r))]
  
  cat("\nB家族内部中位r:", round(median(b_only[-1]), 3), "\n")
  cat("参照基因中位r:", round(median(ref_only), 3), "\n")
  
  if(median(b_only[-1]) < 0.3) {
    cat("🔍 NDUFB7是B家族边缘成员（低连通性）\n")
  } else if(median(b_only[-1]) < 0.5) {
    cat("🔍 NDUFB7是B家族弱-中等成员\n")
  } else {
    cat("🔍 NDUFB7是B家族核心成员\n")
  }
  
  # 与核心/参照比较
  cat("\n与核心亚基比较:\n")
  if("NDUFS4" %in% names(ndufb7_r)) cat("  NDUFS4 (核心): r=", round(ndufb7_r["NDUFS4"], 3), "\n")
  if("NDUFA9" %in% names(ndufb7_r)) cat("  NDUFA9 (Q锚点): r=", round(ndufb7_r["NDUFA9"], 3), "\n")
  if("NDUFV1" %in% names(ndufb7_r)) cat("  NDUFV1 (黄素蛋白): r=", round(ndufb7_r["NDUFV1"], 3), "\n")
  if("NDUFAF2" %in% names(ndufb7_r)) cat("  NDUFAF2 (组装因子): r=", round(ndufb7_r["NDUFAF2"], 3), "\n")
}

# ---------- 7. 生成热图数据 ----------
cor_df <- melt(cor_mat)
names(cor_df) <- c("Gene1", "Gene2", "Correlation")
cor_df$Gene1 <- factor(cor_df$Gene1, levels=rownames(cor_mat))
cor_df$Gene2 <- factor(cor_df$Gene2, levels=rownames(cor_mat))

# 标记B家族
cor_df$Type <- ifelse(grepl("^NDUFB[0-9]+$", cor_df$Gene1) & grepl("^NDUFB[0-9]+$", cor_df$Gene2),
                      "B-B", ifelse(grepl("^NDUFB[0-9]+$", cor_df$Gene1) | grepl("^NDUFB[0-9]+$", cor_df$Gene2),
                                    "B-Ref", "Ref-Ref"))

fwrite(cor_df, "~/Projects/NDUFB7_HF_2026_04_20/03_results/nduf_family_correlation_long.txt", sep="\t")
cat("\n✅ 长格式相关性已保存（供ggplot2热图使用）\n")

cat("\n完成\n")
