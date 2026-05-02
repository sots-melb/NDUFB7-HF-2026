library(data.table)

cat("=" ,rep("=", 69), "\n", sep="")
cat("PROGENy V3: 全列搜索 + 简化通路分析\n")
cat("=" ,rep("=", 69), "\n", sep="")

# 读取平台文件
plat <- fread("~/Projects/NDUFB7_HF_2026_04_20/08_admin/docs/GPL11532-32230.txt", skip=16, header=TRUE)

# 读取表达
gse57338 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt")
ndufb7 <- fread("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE57338/NDUFB7_expression_series_matrix.txt")
ndufb7_vals <- as.numeric(unlist(ndufb7[, -1]))

q25 <- quantile(ndufb7_vals, 0.25)
q75 <- quantile(ndufb7_vals, 0.75)
low_mask <- ndufb7_vals <= q25
high_mask <- ndufb7_vals >= q75
cat("低表达:", sum(low_mask), "高表达:", sum(high_mask), "\n")

# 通路基因
pathways <- list(
  OXPHOS = c("NDUFS1","NDUFS2","NDUFV1","NDUFV2","SDHA","UQCRC1","COX4I1","ATP5A1","ATP5B","ATP5C1"),
  Glycolysis = c("HK2","PFKM","PKM","LDHA","SLC2A1","PGK1","ENO1","GAPDH"),
  ROS = c("SOD1","SOD2","CAT","GPX1","GPX4","PRDX3","TXN2","NQO1"),
  Apoptosis = c("BAX","BAK1","CASP9","CASP3","PARP1","CYCS","BCL2","BCL2L1"),
  Hypoxia = c("VEGFA","EGLN1","HIF1A","CA9","PGK1","LDHA","BNIP3"),
  Fibrosis = c("TGFB1","COL1A1","COL1A2","ACTA2","VIM","FN1","POSTN")
)

# 全列搜索找探针
cat("\n【通路分析】\n")
pathway_results <- data.table(Pathway=character(), N_Genes=integer(), Low_Mean=numeric(), High_Mean=numeric(), Pval=numeric())

for(pw_name in names(pathways)) {
  genes <- pathways[[pw_name]]
  found_probes <- c()
  
  for(g in genes) {
    for(col in names(plat)) {
      vals <- as.character(plat[[col]])
      hits <- grep(paste0(" // ", g, " // "), vals, ignore.case=TRUE)
      if(length(hits) > 0) {
        found_probes <- c(found_probes, as.character(plat[[1]][hits[1]]))
        break
      }
    }
    # 如果严格模式失败，尝试宽松
    if(!g %in% found_probes) {
      for(col in names(plat)) {
        vals <- as.character(plat[[col]])
        hits <- grep(g, vals, ignore.case=TRUE)
        if(length(hits) > 0) {
          found_probes <- c(found_probes, as.character(plat[[1]][hits[1]]))
          break
        }
      }
    }
  }
  
  found_probes <- unique(found_probes)
  probe_ids <- as.character(gse57338[[1]])
  matched <- probe_ids %in% found_probes
  
  if(sum(matched) > 0) {
    expr <- as.matrix(gse57338[matched, -1])
    scores <- colMeans(expr, na.rm=TRUE)
    
    low_mean <- mean(scores[low_mask], na.rm=TRUE)
    high_mean <- mean(scores[high_mask], na.rm=TRUE)
    pval <- wilcox.test(scores[low_mask], scores[high_mask])$p.value
    
    pathway_results <- rbind(pathway_results, data.table(
      Pathway = pw_name,
      N_Genes = sum(matched),
      Low_Mean = low_mean,
      High_Mean = high_mean,
      Pval = pval
    ))
    
    cat(pw_name, "(", sum(matched), "genes): Low=", round(low_mean, 3), 
        " High=", round(high_mean, 3), " p=", format(pval, digits=3), "\n")
  } else {
    cat(pw_name, ": 0 genes found\n")
  }
}

fwrite(pathway_results, "~/Projects/NDUFB7_HF_2026_04_20/03_results/virtual_ko_pathway_v3.txt", sep="\t")

cat("\n【结论】\n")
sig <- pathway_results[Pval < 0.05]
if(nrow(sig) > 0) {
  for(i in 1:nrow(sig)) {
    dir <- ifelse(sig$Low_Mean[i] < sig$High_Mean[i], "↓", "↑")
    cat("  ", sig$Pathway[i], dir, "in NDUFB7-low (p=", format(sig$Pval[i], digits=3), ")\n")
  }
} else {
  cat("  无显著差异通路\n")
}

cat("\n✅ 已保存\n")
cat("\n完成\n")
