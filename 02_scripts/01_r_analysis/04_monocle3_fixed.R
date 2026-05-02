#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     Monocle3拟时序分析                                   ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(ggplot2)
})

srt <- readRDS("03_results/09_single_cell/12_srt_with_stemness.rds")
cat("✅ 加载对象:", ncol(srt), "cells\n")

cat("\n========== 1. Seurat → Monocle3转换 ==========\n")
cds <- as.cell_data_set(srt)
cat("✅ CellDataSet创建:", ncol(cds), "cells\n")

cat("\n========== 2. 预处理与降维 ==========\n")
cds <- preprocess_cds(cds, num_dim = 20)
cat("✅ 预处理完成\n")

cds <- align_cds(cds, alignment_group = "condition")
cat("✅ 批次校正完成\n")

cds <- reduce_dimension(cds)
cat("✅ UMAP降维完成\n")

cds <- cluster_cells(cds)
cat("✅ 聚类完成\n")

cat("\n========== 3. 学习轨迹 ==========\n")
cds <- learn_graph(cds)
cat("✅ 轨迹图学习完成\n")

cat("\n========== 4. 设置根节点 ==========\n")
ctrl_cells <- colnames(cds)[cds$condition == "Control"]
if (length(ctrl_cells) > 0) {
  tryCatch({
    cds <- order_cells(cds, root_cells = ctrl_cells[1:min(10, length(ctrl_cells))])
    cat("✅ 根节点设置完成（Control细胞为起点）\n")
  }, error = function(e) {
    cat("⚠️ root_cells设置失败，自动选择:", conditionMessage(e), "\n")
    cds <<- order_cells(cds)
  })
} else {
  cds <- order_cells(cds)
  cat("⚠️ 无Control细胞，自动选择根节点\n")
}

cat("\n========== 5. Pseudotime统计 ==========\n")
pt <- pseudotime(cds)
cat("Pseudotime范围:", paste(round(range(pt, na.rm=TRUE), 2), collapse=" ~ "), "\n")

# 回写Seurat
srt$pseudotime <- NA
common <- intersect(colnames(srt), colnames(cds))
srt$pseudotime[match(common, colnames(srt))] <- pt[common]

cat("\n========== 6. NDUFB7沿轨迹分析 ==========\n")
if ("NDUFB7" %in% rownames(srt)) {
  ndufb7_expr <- FetchData(srt, vars = "NDUFB7")[,1]
  pt_df <- data.frame(
    pseudotime = srt$pseudotime,
    NDUFB7 = ndufb7_expr,
    Condition = srt$condition,
    Stemness = srt$stemness_score,
    Cluster = srt$seurat_clusters
  )
  pt_df <- pt_df[!is.na(pt_df$pseudotime), ]
  
  ct <- cor.test(pt_df$pseudotime, pt_df$NDUFB7, method = "spearman")
  cat("NDUFB7 vs Pseudotime: rho =", round(ct$estimate, 4), ", p =", format(ct$p.value, digits=4), "\n")
  
  if (ct$estimate < -0.1 && ct$p.value < 0.1) {
    cat("✅ NDUFB7随拟时序下降（支持：远离Control，NDUFB7降低 = 去分化方向）\n")
  } else if (ct$estimate > 0.1) {
    cat("⚠️ NDUFB7随拟时序上升\n")
  } else {
    cat("⚠️ 无显著趋势\n")
  }
  
  cat("\n========== 7. 可视化 ==========\n")
  
  # NDUFB7轨迹
  p1 <- ggplot(pt_df, aes(x = pseudotime, y = NDUFB7, color = Condition)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "loess", se = TRUE, aes(group = 1), color = "black", linetype = "dashed") +
    labs(title = "NDUFB7 along Pseudotime",
         subtitle = paste0("Spearman rho = ", round(ct$estimate, 3)),
         x = "Pseudotime", y = "NDUFB7 Expression") +
    theme_bw() +
    scale_color_manual(values = c("Control" = "#4ECDC4", "HF" = "#FF6B6B"))
  ggsave("03_results/09_single_cell/17_NDUFB7_pseudotime.pdf", p1, width = 8, height = 6)
  ggsave("03_results/09_single_cell/17_NDUFB7_pseudotime.png", p1, width = 8, height = 6, dpi = 300)
  cat("[保存] 17_NDUFB7_pseudotime.pdf/png\n")
  
  # Stemness轨迹（如果有）
  if (sum(!is.na(pt_df$Stemness)) > 10) {
    ct2 <- cor.test(pt_df$pseudotime, pt_df$Stemness, method = "spearman")
    cat("Stemness vs Pseudotime: rho =", round(ct2$estimate, 4), ", p =", format(ct2$p.value, digits=4), "\n")
    
    p2 <- ggplot(pt_df, aes(x = pseudotime, y = Stemness, color = Condition)) +
      geom_point(alpha = 0.3, size = 0.5) +
      geom_smooth(method = "loess", se = TRUE, aes(group = 1), color = "black", linetype = "dashed") +
      labs(title = "Stemness along Pseudotime",
           subtitle = paste0("Spearman rho = ", round(ct2$estimate, 3)),
           x = "Pseudotime", y = "Stemness Score") +
      theme_bw() +
      scale_color_manual(values = c("Control" = "#4ECDC4", "HF" = "#FF6B6B"))
    ggsave("03_results/09_single_cell/18_stemness_pseudotime.pdf", p2, width = 8, height = 6)
    ggsave("03_results/09_single_cell/18_stemness_pseudotime.png", p2, width = 8, height = 6, dpi = 300)
    cat("[保存] 18_stemness_pseudotime.pdf/png\n")
  }
  
  # Monocle3原生图
  pdf("03_results/09_single_cell/19_monocle3_trajectory.pdf", width = 10, height = 8)
  print(plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE) +
    ggtitle("Pseudotime Trajectory"))
  dev.off()
  cat("[保存] 19_monocle3_trajectory.pdf\n")
  
  pdf("03_results/09_single_cell/20_monocle3_condition.pdf", width = 10, height = 8)
  print(plot_cells(cds, color_cells_by = "condition", label_groups_by_cluster = FALSE) +
    ggtitle("Trajectory by Condition"))
  dev.off()
  cat("[保存] 20_monocle3_condition.pdf\n")
  
  # NDUFB7在轨迹上的Monocle3原生图
  if ("NDUFB7" %in% rownames(cds)) {
    pdf("03_results/09_single_cell/20b_monocle3_ndufb7.pdf", width = 10, height = 8)
    print(plot_cells(cds, genes = "NDUFB7", label_groups_by_cluster = FALSE) +
      ggtitle("NDUFB7 on Trajectory"))
    dev.off()
    cat("[保存] 20b_monocle3_ndufb7.pdf\n")
  }
  
} else {
  cat("❌ NDUFB7不在数据中\n")
}

cat("\n========== 8. 保存 ==========\n")
saveRDS(cds, "03_results/09_single_cell/21_monocle3_cds.rds")
cat("[保存] 21_monocle3_cds.rds\n")
saveRDS(srt, "03_results/09_single_cell/22_srt_with_pseudotime.rds")
cat("[保存] 22_srt_with_pseudotime.rds\n")

cat("\n🎉 Monocle3拟时序分析完成！\n")
