#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)
  library(survminer)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V136_GSE57338_Bulk")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V136: GSE57338 Gene-Level铁死亡评分 + Cox")
message("========================================")

# --- 加载GSE57338 gene-level矩阵 ---
rds_file <- "03_results/V133_GSE57338/GSE57338_gene_level.rds"
if (!file.exists(rds_file)) {
  rds_file <- list.files(c("01_data","03_results"), pattern="GSE57338.*\\.rds$", full.names=TRUE, recursive=TRUE)[1]
}
obj <- readRDS(rds_file)
exprs <- obj$exprs
message("[PASS] GSE57338 exprs: ", nrow(exprs), " genes × ", ncol(exprs), " samples")
message("[PASS] NDUFB7 mean: ", round(mean(exprs["NDUFB7", ]), 2))

# --- [1/4] 提取铁死亡基因 ---
ferro_genes <- c("FTL", "SAT1", "NFE2L2", "SLC7A11", "ACSL4", "GPX4", "FTH1")
avail_ferro <- intersect(ferro_genes, rownames(exprs))
message("\n[1/4] 铁死亡基因可用: ", length(avail_ferro), "/", length(ferro_genes))
message("  ", paste(avail_ferro, collapse=", "))

if (length(avail_ferro) >= 3) {
  ferro_mat <- exprs[avail_ferro, ]
  # 标准化评分（z-score per gene, then mean）
  ferro_score <- colMeans(scale(t(ferro_mat)), na.rm = TRUE)
  
  # NDUFB7分层
  ndufb7_expr <- as.numeric(exprs["NDUFB7", ])
  ndufb7_group <- ifelse(ndufb7_expr < median(ndufb7_expr), "NDUFB7_Low", "NDUFB7_High")
  
  # 相关性
  cor_test <- cor.test(ndufb7_expr, ferro_score, method = "spearman")
  message("\n  NDUFB7 vs Ferroptosis Score: rho=", round(cor_test$estimate, 3), 
          " p=", format(cor_test$p.value, digits=2, scientific=TRUE))
  
  # 分组比较
  low_score <- ferro_score[ndufb7_group == "NDUFB7_Low"]
  high_score <- ferro_score[ndufb7_group == "NDUFB7_High"]
  tt <- t.test(low_score, high_score)
  message("  Low vs High NDUFB7: t=", round(tt$statistic, 2), " p=", format(tt$p.value, digits=2, scientific=TRUE))
  
  # 保存
  df_bulk <- data.frame(
    sample = colnames(exprs),
    NDUFB7 = ndufb7_expr,
    ferroptosis_score = ferro_score,
    NDUFB7_group = ndufb7_group
  )
  write.csv(df_bulk, file.path(outdir, "V136_GSE57338_bulk_ferroptosis.csv"), row.names=FALSE)
  
  # 可视化
  p <- ggplot(df_bulk, aes(x=NDUFB7_group, y=ferroptosis_score, fill=NDUFB7_group)) +
    geom_boxplot(alpha=0.7) + geom_jitter(width=0.2, alpha=0.3) +
    scale_fill_manual(values=c("NDUFB7_Low"="#440154", "NDUFB7_High"="#35B779")) +
    labs(title="Bulk Ferroptosis Score by NDUFB7 Level (GSE57338)", 
         subtitle=paste0("Spearman rho=", round(cor_test$estimate,3), " p=", format(cor_test$p.value, digits=1, scientific=TRUE)),
         x="NDUFB7 Group", y="Ferroptosis Score (z-mean)") +
    theme_minimal()
  ggsave(file.path(outdir, "V136_bulk_ferroptosis_boxplot.png"), p, width=6, height=5, dpi=300)
  
  message("[PASS] Bulk ferroptosis scoring complete")
} else {
  message("[WARN] Insufficient ferroptosis genes available")
}

# --- [2/4] ACSL4/GPX4比值（Bulk层面）---
if (all(c("ACSL4", "GPX4") %in% rownames(exprs))) {
  ratio_bulk <- as.numeric(exprs["ACSL4", ]) / (as.numeric(exprs["GPX4", ]) + 1e-6)
  cor_ratio <- cor.test(ndufb7_expr, ratio_bulk, method = "spearman")
  message("\n[2/4] ACSL4/GPX4 ratio vs NDUFB7: rho=", round(cor_ratio$estimate, 3), 
          " p=", format(cor_ratio$p.value, digits=2, scientific=TRUE))
  
  df_bulk$ACSL4_GPX4_ratio <- ratio_bulk
  write.csv(df_bulk, file.path(outdir, "V136_GSE57338_bulk_ferroptosis.csv"), row.names=FALSE)
} else {
  message("[SKIP] ACSL4/GPX4 not both available")
}

# --- [3/4] Cox生存分析（如果有临床数据）---
clinical_candidates <- list.files(c("01_data","03_results"), pattern="GSE57338.*clinical|GSE57338.*pheno|57338.*survival", full.names=TRUE, recursive=TRUE)
clinical_file <- clinical_candidates[1]

if (!is.na(clinical_file) && file.exists(clinical_file)) {
  message("\n[3/4] Clinical data found: ", basename(clinical_file))
  # 简化：假设有time/status列，或需根据实际格式调整
  message("[INFO] Cox analysis requires time/status columns — please verify clinical file format")
  # 框架代码，实际执行需匹配列名
} else {
  message("\n[SKIP] No clinical data found for Cox analysis")
  message("[ACTION] GSE57338 clinical data may be in series matrix or separate file")
}

# --- [4/4] 保存 ---
message("\n[4/4] 保存汇总")
write.csv(data.frame(
  Metric = c("N_samples", "N_genes", "NDUFB7_mean", "NDUFB7_median", "NDUFB7_sd",
             "Ferro_genes_used", "Spearman_rho", "Spearman_p", "T_test_p"),
  Value = c(ncol(exprs), nrow(exprs), round(mean(ndufb7_expr),2), round(median(ndufb7_expr),2), round(sd(ndufb7_expr),2),
            length(avail_ferro), round(cor_test$estimate,3), format(cor_test$p.value, scientific=TRUE), format(tt$p.value, scientific=TRUE))
), file.path(outdir, "V136_summary_stats.csv"), row.names=FALSE)

message("[DONE] V136: ", outdir)
