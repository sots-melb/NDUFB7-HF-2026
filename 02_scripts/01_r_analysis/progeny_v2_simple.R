library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("PROGENy简化版: 直接使用GSE57338差异基因\n")
cat("=" ,rep("=", 69), "\n", sep="")

# ---------- 1. 读取NDUFB7和表型 ----------
ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
ndufb7_vals <- as.numeric(unlist(ndufb7[, -1]))
pheno <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_phenotype.txt")

# ---------- 2. 定义NDUFB7低表达组 ----------
q25 <- quantile(ndufb7_vals, 0.25)
q75 <- quantile(ndufb7_vals, 0.75)
low_mask <- ndufb7_vals <= q25
high_mask <- ndufb7_vals >= q75

cat("低表达:", sum(low_mask), "高表达:", sum(high_mask), "\n")

# ---------- 3. 手动通路分析（不依赖PROGENy包）----------
# 使用已知的NDUFB7相关通路基因
pathways <- list(
  OXPHOS = c("NDUFS1","NDUFS2","NDUFV1","NDUFV2","SDHA","UQCRC1","COX4I1","ATP5A1","ATP5B"),
  Glycolysis = c("HK2","PFKM","PKM","LDHA","SLC2A1","PGK1","ENO1"),
  ROS_Defense = c("SOD1","SOD2","CAT","GPX1","GPX4","PRDX3","TXN2"),
  Apoptosis = c("BAX","BAK1","CASP9","CASP3","PARP1","CYCS","BCL2"),
  Mitophagy = c("PINK1","PRKN","BNIP3","BNIP3L","SQSTM1"),
  Hypoxia = c("VEGFA","EGLN1","HIF1A","CA9","PGK1","LDHA"),
  Fibrosis = c("TGFB1","COL1A1","COL1A2","ACTA2","VIM","FN1")
)

# 读取表达矩阵（使用探针ID，不映射gene symbol）
gse57338 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")
sample_cols <- setdiff(names(gse57338), names(gse57338)[1])
expr_mat <- as.matrix(gse57338[, ..sample_cols])

# 读取平台文件建立映射（简单版：只找通路基因）
plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)
anno_col <- names(plat)[ncol(plat)]

# 为每个通路找探针
cat("\n【通路分析】\n")
pathway_results <- data.table(Pathway=character(), Low_Mean=numeric(), High_Mean=numeric(), Pval=numeric(), N_Genes=integer())

for(pw_name in names(pathways)) {
  genes <- pathways[[pw_name]]
  
  # 在平台注释中搜索这些基因
  found_probes <- c()
  for(g in genes) {
    pattern <- paste0(" // ", g, " //")
    hits <- plat[grepl(pattern, plat[[anno_col]], ignore.case=TRUE)]
    if(nrow(hits) > 0) {
      found_probes <- c(found_probes, as.character(hits[[1]][1]))  # 取第一个探针
    }
  }
  found_probes <- unique(found_probes)
  
  if(length(found_probes) > 0) {
    # 在表达矩阵中提取
    probe_mask <- rownames(expr_mat) %in% found_probes
    if(sum(probe_mask) > 0) {
      scores <- colMeans(expr_mat[probe_mask, , drop=FALSE], na.rm=TRUE)
      
      low_mean <- mean(scores[low_mask], na.rm=TRUE)
      high_mean <- mean(scores[high_mask], na.rm=TRUE)
      pval <- wilcox.test(scores[low_mask], scores[high_mask])$p.value
      
      pathway_results <- rbind(pathway_results, data.table(
        Pathway = pw_name,
        Low_Mean = low_mean,
        High_Mean = high_mean,
        Pval = pval,
        N_Genes = sum(probe_mask)
      ))
      
      cat(pw_name, "(", sum(probe_mask), "genes): Low=", round(low_mean, 3), 
          " High=", round(high_mean, 3), " p=", format(pval, digits=3), "\n")
    }
  }
}

fwrite(pathway_results, "~/Projects/NDUFB7_HF_2026_04_20/03_results/virtual_ko_pathway_v2.txt", sep="\t")

cat("\n【虚拟KO结论】\n")
sig <- pathway_results[Pval < 0.05]
if(nrow(sig) > 0) {
  cat("显著差异通路 (p<0.05):\n")
  for(i in 1:nrow(sig)) {
    direction <- ifelse(sig$Low_Mean[i] < sig$High_Mean[i], "↓", "↑")
    cat("  ", sig$Pathway[i], direction, "in NDUFB7-low (p=", format(sig$Pval[i], digits=3), ")\n")
  }
} else {
  cat("  无显著差异通路（样本量或基因覆盖不足）\n")
}

cat("\n✅ 结果已保存\n")
cat("\n完成\n")
