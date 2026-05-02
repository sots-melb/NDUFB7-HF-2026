#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)
cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     StemID 细胞干性推断分析 (Seurat V4兼容版)            ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

cat("\n========== 0. 加载Seurat对象 ==========\n")
srt <- readRDS("03_results/09_single_cell/06_srt_v4_processed.rds")
cat("✅ 加载对象:", ncol(srt), "cells x", nrow(srt), "genes\n")
cat("Condition分布:\n")
print(table(srt$condition))

cat("\n========== 1. Cluster 4关键发现预览 ==========\n")
ndufb7_cluster <- aggregate(NDUFB7 ~ seurat_clusters + condition,
                             data = FetchData(srt, vars = c("NDUFB7","seurat_clusters","condition")),
                             FUN = mean)
cat("各簇NDUFB7平均表达:\n")
print(ndufb7_cluster[order(-ndufb7_cluster$NDUFB7), ])
cat("\n🎯 Cluster 4: HF特异, NDUFB7最低\n")

cat("\n========== 2. 加载RaceID ==========\n")
if (!requireNamespace("RaceID", quietly = TRUE)) {
  cat("⏳ RaceID未安装，正在安装（约2-5分钟）...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  BiocManager::install("RaceID", update = FALSE, ask = FALSE, quiet = TRUE)
}
library(RaceID)
cat("✅ RaceID加载成功\n")

cat("\n========== 3. 准备RaceID输入 ==========\n")
expr <- GetAssayData(srt, assay = "RNA", slot = "counts")
cat("Counts矩阵:", nrow(expr), "x", ncol(expr), "\n")
expr_dense <- as.matrix(expr)
cat("✅ 转换为dense matrix\n")

cat("\n========== 4. StemID计算（约2-5分钟）==========\n")
sc <- SCseq(expr_dense)
cat("SCseq对象创建完成\n")
sc <- filterdata(sc, mintotal = 3000, minexpr = 5, minnumber = 5,
                 maxexpr = Inf, downsample = FALSE)
cat("过滤后细胞数:", ncol(sc@ndata), "\n")
ltr <- Ltree(sc)
cat("谱系树创建完成\n")
cat("计算熵（干性指数）...\n")
ltr <- compentropy(ltr)
stemness_scores <- ltr@entropy
cat("✅ 干性计算完成\n")
cat("Stemness范围:", paste(round(range(stemness_scores, na.rm=TRUE), 4), collapse=" ~ "), "\n")

cat("\n========== 5. 结果回写Seurat ==========\n")
common_cells <- intersect(colnames(srt), names(stemness_scores))
cat("匹配细胞:", length(common_cells), "/", ncol(srt), "\n")
srt$stemness_score <- NA
srt$stemness_score[match(common_cells, colnames(srt))] <- stemness_scores[common_cells]

stem_df <- srt@meta.data %>%
  group_by(seurat_clusters, condition) %>%
  summarise(
    stemness_mean = mean(stemness_score, na.rm = TRUE),
    stemness_median = median(stemness_score, na.rm = TRUE),
    n = sum(!is.na(stemness_score)),
    .groups = "drop"
  )
cat("\n各簇Stemness统计:\n")
print(as.data.frame(stem_df))

cat("\n========== 6. NDUFB7与干性关联分析 ==========\n")
ndufb7_expr <- FetchData(srt, vars = "NDUFB7")[,1]
valid <- !is.na(srt$stemness_score) & !is.na(ndufb7_expr)
cat("有效样本:", sum(valid), "\n")

cor_result <- NULL
if (sum(valid) > 10) {
  cor_result <- cor.test(srt$stemness_score[valid], ndufb7_expr[valid], method = "spearman")
  cat("\n✅✅✅ Spearman相关结果 ✅✅✅\n")
  cat("rho =", round(cor_result$estimate, 4), "\n")
  cat("p =", format(cor_result$p.value, digits = 4), "\n")
  cat("n =", sum(valid), "\n")
  
  if (cor_result$estimate < -0.2 && cor_result$p.value < 0.05) {
    cat("\n🎉🎉🎉 理想结果：NDUFB7低表达与干性升高显著负相关！\n")
    cat("强力支持'NDUFB7下调→线粒体功能障碍→心肌细胞去分化→心衰'因果链\n")
  } else if (cor_result$estimate < -0.1 && cor_result$p.value < 0.1) {
    cat("\n✅ 趋势性证据：负相关趋势明显\n")
  } else if (cor_result$estimate > 0) {
    cat("\n⚠️ 意外正相关，需调整叙事\n")
  } else {
    cat("\n⚠️ 结果不显著\n")
  }
  
  stats <- data.frame(
    rho = round(cor_result$estimate, 4),
    p_value = format(cor_result$p.value, digits = 4),
    n_cells = sum(valid),
    direction = ifelse(cor_result$estimate < 0, "Negative", "Positive"),
    interpretation = ifelse(cor_result$estimate < -0.2 && cor_result$p.value < 0.05,
                            "Strong_support",
                            ifelse(cor_result$estimate < -0.1 && cor_result$p.value < 0.1,
                                   "Trend_support", "Weak_or_no_support"))
  )
  write.csv(stats, "03_results/09_single_cell/08_stemid_correlation.csv", row.names = FALSE)
  cat("[保存] 08_stemid_correlation.csv\n")
  
  c4_cells <- srt$seurat_clusters == 4
  if (sum(c4_cells, na.rm=TRUE) > 5) {
    cat("\nCluster 4特异性:\n")
    cat("  Stemness均值:", round(mean(srt$stemness_score[c4_cells], na.rm=TRUE), 4), "\n")
    cat("  NDUFB7均值:", round(mean(ndufb7_expr[c4_cells], na.rm=TRUE), 4), "\n")
    cat("  细胞数:", sum(c4_cells), "(HF:", sum(c4_cells & srt$condition=="HF"), ")\n")
  }
} else {
  cat("❌ 有效样本不足\n")
}

cat("\n========== 7. 可视化 ==========\n")
plot_data <- data.frame(
  NDUFB7 = ndufb7_expr,
  Stemness = srt$stemness_score,
  Condition = srt$condition,
  Cluster = srt$seurat_clusters,
  stringsAsFactors = FALSE
)

p1 <- ggplot(plot_data, aes(x = NDUFB7, y = Stemness, color = Condition)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, aes(group = 1), color = "black", linetype = "dashed") +
  labs(title = "NDUFB7 vs Stemness Score",
       subtitle = paste0("Spearman rho = ", round(cor_result$estimate, 3),
                       ", p = ", format(cor_result$p.value, digits = 3)),
       x = "NDUFB7 Expression", y = "Stemness Score") +
  theme_bw() +
  scale_color_manual(values = c("Control" = "#4ECDC4", "HF" = "#FF6B6B"))
ggsave("03_results/09_single_cell/09_NDUFB7_stemness_scatter.pdf", p1, width = 8, height = 6)
cat("[保存] 09_NDUFB7_stemness_scatter.pdf\n")

pdf("03_results/09_single_cell/10_UMAP_stemness.pdf", width = 8, height = 6)
print(FeaturePlot(srt, features = "stemness_score", pt.size = 0.5) +
  ggtitle("Stemness Score on UMAP") +
  scale_color_gradientn(colors = c("blue", "green", "yellow", "red")))
dev.off()
cat("[保存] 10_UMAP_stemness.pdf\n")

pdf("03_results/09_single_cell/11_stemness_by_cluster.pdf", width = 10, height = 6)
print(VlnPlot(srt, features = "stemness_score", pt.size = 0.1, group.by = "seurat_clusters") +
  ggtitle("Stemness Score by Cluster"))
dev.off()
cat("[保存] 11_stemness_by_cluster.pdf\n")

saveRDS(srt, "03_results/09_single_cell/12_srt_with_stemness.rds")
cat("\n[保存] 12_srt_with_stemness.rds\n")
cat("\n🎉 StemID分析完成！\n\n")
