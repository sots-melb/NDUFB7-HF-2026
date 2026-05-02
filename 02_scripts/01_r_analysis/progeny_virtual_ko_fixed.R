library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("PROGENy虚拟敲除（修复版）: 探针ID→Gene Symbol映射\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 1. 读取平台文件建立探针→基因映射 ----------
plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)
anno_col <- names(plat)[ncol(plat)]

get_gene <- function(x) {
  m <- regmatches(x, gregexpr("// [A-Z0-9]+ //", x))
  if(length(m[[1]]) > 0) gsub(" // ", "", m[[1]][1]) else NA
}

plat$GeneSymbol <- sapply(plat[[anno_col]], get_gene)
plat <- plat[!is.na(plat$GeneSymbol)]
probe_to_gene <- setNames(plat$GeneSymbol, as.character(plat[[1]]))
cat("探针→基因映射:", length(probe_to_gene), "\n")

# ---------- 2. 读取GSE57338并映射基因名 ----------
gse57338 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")
sample_cols <- setdiff(names(gse57338), names(gse57338)[1])

expr_mat <- as.matrix(gse57338[, ..sample_cols])
rownames(expr_mat) <- as.character(gse57338[[1]])

# 映射到gene symbol
gene_symbols <- probe_to_gene[rownames(expr_mat)]
# 去重：同一基因多个探针取均值
unique_genes <- unique(na.omit(gene_symbols))
cat("唯一基因:", length(unique_genes), "\n")

gene_expr <- sapply(unique_genes, function(g) {
  probes <- names(gene_symbols)[gene_symbols == g]
  if(length(probes) > 0) {
    colMeans(expr_mat[probes, , drop=FALSE], na.rm=TRUE)
  } else rep(NA, ncol(expr_mat))
})
gene_expr <- t(gene_expr)
colnames(gene_expr) <- sample_cols
cat("Gene Symbol表达矩阵:", nrow(gene_expr), "genes x", ncol(gene_expr), "samples\n")

# ---------- 3. NDUFB7高低分组 ----------
ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
ndufb7_vals <- as.numeric(unlist(ndufb7[, -1]))

q25 <- quantile(ndufb7_vals, 0.25)
q75 <- quantile(ndufb7_vals, 0.75)
low_mask <- ndufb7_vals <= q25
high_mask <- ndufb7_vals >= q75

cat("\n低表达(bottom 25%):", sum(low_mask), " 高表达(top 25%):", sum(high_mask), "\n")

# ---------- 4. PROGENy分析（使用gene symbol矩阵）----------
if(require("progeny", quietly=TRUE)) {
  cat("\n【PROGENy通路预测】\n")
  
  # 子集：低+高表达样本
  sub_expr <- gene_expr[, low_mask | high_mask]
  # 去除NA行
  sub_expr <- sub_expr[complete.cases(sub_expr), ]
  cat("PROGENy输入:", nrow(sub_expr), "genes x", ncol(sub_expr), "samples\n")
  
  progeny_scores <- progeny(sub_expr, scale=TRUE, organism="Human", top=500)
  cat("\n通路活性评分（均值）:\n")
  print(round(colMeans(progeny_scores), 3))
  
  # 比较低vs高
  group <- ifelse(colMeans(sub_expr)["NDUFB7"] <= q25, "Low", "High")  # 简化分组
  # 实际上应该用原始分组
  low_scores <- progeny_scores[, low_mask[low_mask | high_mask]]
  high_scores <- progeny_scores[, high_mask[low_mask | high_mask]]
  
  cat("\n低表达组通路活性:\n")
  print(round(rowMeans(low_scores), 3))
  cat("\n高表达组通路活性:\n")
  print(round(rowMeans(high_scores), 3))
  
  # 差异检验
  cat("\n【通路差异】\n")
  for(pw in rownames(progeny_scores)) {
    wt <- wilcox.test(low_scores[pw, ], high_scores[pw, ])
    cat(pw, ": Low=", round(mean(low_scores[pw, ]), 3), 
        " High=", round(mean(high_scores[pw, ]), 3),
        " p=", format(wt$p.value, digits=3), "\n")
  }
  
  fwrite(as.data.table(progeny_scores, keep.rownames="Pathway"),
         "~/Projects/NDUFB7_HF_2026_04_20/03_results/progeny_scores_fixed.txt", sep="\t")
  cat("\n✅ PROGENy结果已保存\n")
  
} else {
  cat("\n⚠️ PROGENy不可用，使用手动通路分析\n")
  
  # 手动定义通路基因
  pathways <- list(
    OXPHOS = c("NDUFS1","NDUFS2","NDUFV1","SDHA","UQCRC1","COX4I1","ATP5A1"),
    Hypoxia = c("VEGFA","EGLN1","PGK1","LDHA","BNIP3","CA9"),
    Apoptosis = c("BAX","BAK1","CASP9","PARP1","CYCS"),
    P53 = c("TP53","CDKN1A","MDM2","GADD45A"),
    TNFa = c("TNF","NFKB1","RELA","ICAM1","VCAM1"),
    ROS = c("SOD1","SOD2","CAT","GPX1","PRDX3")
  )
  
  cat("\n【简化通路分析】\n")
  for(pw_name in names(pathways)) {
    genes <- pathways[[pw_name]]
    pw_mask <- rownames(gene_expr) %in% genes
    if(sum(pw_mask) > 0) {
      scores <- colMeans(gene_expr[pw_mask, ], na.rm=TRUE)
      cat(pw_name, "(n=", sum(pw_mask), "genes): ",
          "Low=", round(mean(scores[low_mask]), 3),
          " High=", round(mean(scores[high_mask]), 3),
          " p=", format(wilcox.test(scores[low_mask], scores[high_mask])$p.value, digits=3), "\n")
    }
  }
}

cat("\n完成\n")
