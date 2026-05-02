#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     NDUFB7-NDUFB8共表达分析 (B亚家族模块)                ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
})

srt <- readRDS("03_results/09_single_cell/12_srt_with_stemness.rds")
cat("✅ 加载对象:", ncol(srt), "cells\n")

b_family <- c("NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11")
available <- b_family[b_family %in% rownames(srt)]
cat("B亚家族可用基因:", paste(available, collapse = ", "), "\n")

if (length(available) >= 2) {
  expr_mat <- as.matrix(GetAssayData(srt, assay = "RNA", slot = "data")[available, ])
  
  cat("\n========== B亚家族Spearman相关矩阵 ==========\n")
  cor_mat <- cor(t(expr_mat), method = "spearman")
  print(round(cor_mat, 3))
  
  # 保存热图
  pdf("03_results/09_single_cell/14_B_family_cor_heatmap.pdf", width = 7, height = 6)
  pheatmap(cor_mat, display_numbers = TRUE, number_format = "%.3f",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "B-subfamily Co-expression (Spearman)")
  dev.off()
  cat("[保存] 14_B_family_cor_heatmap.pdf\n")
  
  png("03_results/09_single_cell/14_B_family_cor_heatmap.png", width = 700, height = 600, res = 100)
  pheatmap(cor_mat, display_numbers = TRUE, number_format = "%.3f",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "B-subfamily Co-expression (Spearman)")
  dev.off()
  cat("[保存] 14_B_family_cor_heatmap.png\n")
  
  # NDUFB7-NDUFB8 特异性检验
  if ("NDUFB7" %in% available && "NDUFB8" %in% available) {
    ndufb7_expr <- as.numeric(expr_mat["NDUFB7", ])
    ndufb8_expr <- as.numeric(expr_mat["NDUFB8", ])
    ct <- cor.test(ndufb7_expr, ndufb8_expr, method = "spearman")
    
    cat("\n✅✅✅ NDUFB7-NDUFB8共表达 ✅✅✅\n")
    cat("Spearman rho =", round(ct$estimate, 4), "\n")
    cat("p =", format(ct$p.value, digits = 4), "\n")
    
    interp <- ifelse(ct$estimate > 0.3 && ct$p.value < 0.05, "Co_module_support",
                     ifelse(ct$estimate > 0.1 && ct$p.value < 0.1, "Trend_support", "Weak"))
    if (ct$estimate > 0.3 && ct$p.value < 0.05) {
      cat("🎉 显著正相关！支持B亚家族'共失调模块'假说\n")
    } else if (ct$estimate > 0.1) {
      cat("✅ 趋势性正相关\n")
    } else {
      cat("⚠️ 相关较弱或负相关\n")
    }
    
    write.csv(data.frame(
      pair = "NDUFB7-NDUFB8",
      rho = round(ct$estimate, 4),
      p_value = format(ct$p.value, digits = 4),
      n_cells = ncol(expr_mat),
      interpretation = interp
    ), "03_results/09_single_cell/13_ndufb7_ndufb8_correlation.csv", row.names = FALSE)
    cat("[保存] 13_ndufb7_ndufb8_correlation.csv\n")
    
    # 散点图
    plot_df <- data.frame(
      NDUFB7 = ndufb7_expr,
      NDUFB8 = ndufb8_expr,
      Condition = srt$condition,
      Cluster = srt$seurat_clusters,
      stringsAsFactors = FALSE
    )
    
    p <- ggplot(plot_df, aes(x = NDUFB7, y = NDUFB8, color = Condition)) +
      geom_point(alpha = 0.3, size = 0.5) +
      geom_smooth(method = "lm", se = TRUE, aes(group = 1), color = "black", linetype = "dashed") +
      labs(title = "NDUFB7 vs NDUFB8 Co-expression",
           subtitle = paste0("Spearman rho = ", round(ct$estimate, 3),
                           ", p = ", format(ct$p.value, digits = 3)),
           x = "NDUFB7 Expression", y = "NDUFB8 Expression") +
      theme_bw() +
      scale_color_manual(values = c("Control" = "#4ECDC4", "HF" = "#FF6B6B"))
    ggsave("03_results/09_single_cell/15_NDUFB7_vs_NDUFB8_scatter.pdf", p, width = 8, height = 6)
    ggsave("03_results/09_single_cell/15_NDUFB7_vs_NDUFB8_scatter.png", p, width = 8, height = 6, dpi = 300)
    cat("[保存] 15_NDUFB7_vs_NDUFB8_scatter.pdf/png\n")
  }
  
  # 按Cluster的平均表达热图
  cat("\n========== B亚家族按Cluster表达 ==========\n")
  avg_expr <- AverageExpression(srt, features = available, group.by = "seurat_clusters")$RNA
  print(round(avg_expr, 3))
  
  pdf("03_results/09_single_cell/16_B_family_cluster_heatmap.pdf", width = 8, height = 5)
  pheatmap(avg_expr, scale = "row",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "B-subfamily by Cluster (Z-score)")
  dev.off()
  cat("[保存] 16_B_family_cluster_heatmap.pdf\n")
  
  png("03_results/09_single_cell/16_B_family_cluster_heatmap.png", width = 800, height = 500, res = 100)
  pheatmap(avg_expr, scale = "row",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "B-subfamily by Cluster (Z-score)")
  dev.off()
  cat("[保存] 16_B_family_cluster_heatmap.png\n")
  
} else {
  cat("❌ B亚家族可用成员不足2个，无法分析\n")
}

cat("\n🎉 共表达分析完成！\n")
