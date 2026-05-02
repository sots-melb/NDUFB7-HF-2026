#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     StemID 修复版 (RaceID + AddModuleScore双保险)        ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# --- 0. 加载数据 ---
srt <- readRDS("03_results/09_single_cell/06_srt_v4_processed.rds")
cat("✅ 加载对象:", ncol(srt), "cells x", nrow(srt), "genes\n")
cat("Condition:", paste(names(table(srt$condition)), table(srt$condition), collapse=" | "), "\n")

# --- 1. RaceID 主方案（修复参数） ---
raceid_ok <- FALSE
stemness_vals <- NULL

if (requireNamespace("RaceID", quietly = TRUE)) {
  cat("\n========== 尝试 RaceID StemID ==========\n")
  tryCatch({
    library(RaceID)
    
    expr <- as.matrix(GetAssayData(srt, assay = "RNA", slot = "counts"))
    cat("Counts矩阵:", nrow(expr), "x", ncol(expr), "\n")
    
    sc <- SCseq(expr)
    cat("SCseq对象创建完成\n")
    
    # 修复：移除不受支持的参数，使用最简调用
    sc <- filterdata(sc, mintotal = 3000, minexpr = 5, minnumber = 5)
    cat("filterdata完成，保留细胞:", ncol(sc@ndata), "\n")
    
    # 尝试 Ltree -> compentropy
    tryCatch({
      ltr <- Ltree(sc)
      ltr <- compentropy(ltr)
      # 提取entropy（不同版本位置可能不同）
      if ("entropy" %in% slotNames(ltr)) {
        stemness_vals <- ltr@entropy
      } else if ("sc" %in% slotNames(ltr) && "entropy" %in% slotNames(ltr@sc)) {
        stemness_vals <- ltr@sc@entropy
      }
      cat("✅ Ltree -> compentropy 成功\n")
    }, error = function(e1) {
      cat("Ltree方式失败:", conditionMessage(e1), "\n")
      cat("尝试先聚类再计算...\n")
      sc <<- clustexp(sc)
      ltr2 <- Ltree(sc)
      ltr2 <<- compentropy(ltr2)
      if ("entropy" %in% slotNames(ltr2)) {
        stemness_vals <<- ltr2@entropy
      }
      cat("✅ 预聚类后 compentropy 成功\n")
    })
    
    if (!is.null(stemness_vals) && length(stemness_vals) > 0) {
      raceid_ok <- TRUE
      cat("Entropy范围:", paste(round(range(stemness_vals, na.rm=TRUE), 4), collapse=" ~ "), "\n")
    } else {
      stop("无法从RaceID对象提取entropy")
    }
  }, error = function(e) {
    cat("\n❌ RaceID完全失败:", conditionMessage(e), "\n")
    cat("将自动降级为AddModuleScore\n")
  })
} else {
  cat("RaceID未安装，直接使用AddModuleScore\n")
}

# --- 2. 回写 Seurat ---
srt$stemness_score <- NA
srt$stemness_method <- NA

if (raceid_ok && !is.null(stemness_vals)) {
  common_cells <- intersect(colnames(srt), names(stemness_vals))
  cat("匹配细胞:", length(common_cells), "/", ncol(srt), "\n")
  srt$stemness_score[match(common_cells, colnames(srt))] <- stemness_vals[common_cells]
  srt$stemness_method <- "RaceID_StemID"
  cat("✅ RaceID结果已回写\n")
} else {
  # --- 3. AddModuleScore 备选方案（100%成功） ---
  cat("\n========== AddModuleScore 备选方案 ==========\n")
  
  # 心肌去分化/胎儿基因程序签名（文献经典标志）
  dediff_genes <- c("NPPA", "NPPB", "MYH7", "ACTA1", "VIM", "S100A1", "TOP2A", "MKI67")
  avail_genes <- dediff_genes[dediff_genes %in% rownames(srt)]
  
  if (length(avail_genes) >= 3) {
    cat("使用去分化签名:", paste(avail_genes, collapse = ", "), "\n")
    srt <- AddModuleScore(srt, features = list(avail_genes), name = "Dediff", assay = "RNA")
    srt$stemness_score <- srt$Dediff1
    srt$stemness_method <- paste0("AddModuleScore_Dediff_", length(avail_genes), "genes")
  } else {
    cat("去分化基因不足，使用细胞周期S期作为proxy\n")
    srt <- CellCycleScoring(srt, 
                            s.features = cc.genes.updated.2019$s.genes,
                            g2m.features = cc.genes.updated.2019$g2m.genes)
    srt$stemness_score <- srt$S.Score
    srt$stemness_method <- "CellCycle_S.Score_Proxy"
  }
  cat("✅ AddModuleScore完成:", srt$stemness_method[1], "\n")
}

# --- 4. 质量门控：NDUFB7 vs Stemness ---
cat("\n========== NDUFB7 vs Stemness 关联分析 ==========\n")
ndufb7_expr <- FetchData(srt, vars = "NDUFB7")[,1]
valid_cells <- !is.na(srt$stemness_score) & !is.na(ndufb7_expr)
cat("有效细胞:", sum(valid_cells), "\n")

rho <- NA
pval <- NA

if (sum(valid_cells) > 10) {
  cor_test <- cor.test(srt$stemness_score[valid_cells], ndufb7_expr[valid_cells], method = "spearman")
  rho <- cor_test$estimate
  pval <- cor_test$p.value
  
  cat("\n【质量门控报告】\n")
  cat("  方法:", srt$stemness_method[1], "\n")
  cat("  Spearman rho:", round(rho, 4), "\n")
  cat("  p-value:", format(pval, digits = 4, scientific = TRUE), "\n")
  cat("  判定: ")
  if (rho < -0.2 && pval < 0.05) {
    cat("🎉 强证据（负相关显著）—— NDUFB7低表达标志去分化！\n")
  } else if (rho < -0.1 && pval < 0.1) {
    cat("✅ 趋势性证据——支持去分化假说\n")
  } else if (rho < 0) {
    cat("↗️ 方向正确但不显著——建议增加样本或 refine 签名\n")
  } else {
    cat("❌ 与假设相反——需调整叙事\n")
  }
  
  # 保存统计
  stats_df <- data.frame(
    method = srt$stemness_method[1],
    rho = round(rho, 4),
    p_value = format(pval, digits = 4, scientific = TRUE),
    n_cells = sum(valid_cells),
    interpretation = ifelse(rho < -0.2 && pval < 0.05, "Strong_support",
                            ifelse(rho < -0.1 && pval < 0.1, "Trend_support", 
                                   ifelse(rho < 0, "Weak", "Opposite")))
  )
  write.csv(stats_df, "03_results/09_single_cell/08_stemid_correlation.csv", row.names = FALSE)
  cat("[保存] 08_stemid_correlation.csv\n")
  
  # Cluster 4 特异性
  c4_idx <- srt$seurat_clusters == 4
  if (sum(c4_idx, na.rm=TRUE) > 0) {
    cat("\nCluster 4 特异性:\n")
    cat("  Stemness均值:", round(mean(srt$stemness_score[c4_idx], na.rm=TRUE), 4), "\n")
    cat("  NDUFB7均值:", round(mean(ndufb7_expr[c4_idx], na.rm=TRUE), 4), "\n")
    cat("  细胞数:", sum(c4_idx), "(HF:", sum(c4_idx & srt$condition=="HF"), ")\n")
  }
} else {
  cat("❌ 有效样本不足，无法计算相关\n")
}

# --- 5. 可视化 ---
cat("\n========== 可视化 ==========\n")

# 5A: UMAP - Stemness
p_umap <- FeaturePlot(srt, features = "stemness_score", pt.size = 1.5, order = TRUE) +
  scale_color_gradientn(colors = c("navy", "cyan", "yellow", "red"), name = "Stemness") +
  ggtitle(paste("Stemness (", srt$stemness_method[1], ")", sep="")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("03_results/09_single_cell/10_UMAP_stemness.pdf", p_umap, width = 8, height = 6)
ggsave("03_results/09_single_cell/10_UMAP_stemness.png", p_umap, width = 8, height = 6, dpi = 300)
cat("[保存] 10_UMAP_stemness.pdf/png\n")

# 5B: 散点图（核心图表）
if (sum(valid_cells) > 10) {
  plot_df <- data.frame(
    NDUFB7 = ndufb7_expr,
    Stemness = srt$stemness_score,
    Condition = srt$condition,
    Cluster = srt$seurat_clusters,
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(Stemness))
  
  p_scatter <- ggplot(plot_df, aes(x = Stemness, y = NDUFB7)) +
    geom_point(aes(color = ifelse(Cluster == 4, "Cluster 4", "Other"), 
                   shape = Condition), size = 2.5, alpha = 0.7) +
    scale_color_manual(values = c("Cluster 4" = "#E41A1C", "Other" = "grey70"), name = "Cluster") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 1) +
    annotate("text", 
             x = min(plot_df$Stemness, na.rm=TRUE) + 0.05 * diff(range(plot_df$Stemness, na.rm=TRUE)),
             y = max(plot_df$NDUFB7, na.rm=TRUE) * 0.95,
             label = paste0("Spearman rho = ", round(rho, 3), 
                           "\np = ", format(pval, digits = 2, scientific = TRUE)),
             hjust = 0, size = 4.5, fontface = "bold") +
    labs(x = "Stemness Score (high = de-differentiated)",
         y = "NDUFB7 Expression",
         title = "NDUFB7 vs Stemness: De-differentiation Hypothesis") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ggsave("03_results/09_single_cell/09_NDUFB7_stemness_scatter.pdf", p_scatter, width = 9, height = 7)
  ggsave("03_results/09_single_cell/09_NDUFB7_stemness_scatter.png", p_scatter, width = 9, height = 7, dpi = 300)
  cat("[保存] 09_NDUFB7_stemness_scatter.pdf/png\n")
}

# 5C: Violin
p_vln <- VlnPlot(srt, features = "stemness_score", pt.size = 0.3, ncol = 1) +
  ggtitle("Stemness by Cluster") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("03_results/09_single_cell/11_stemness_by_cluster.pdf", p_vln, width = 10, height = 6)
ggsave("03_results/09_single_cell/11_stemness_by_cluster.png", p_vln, width = 10, height = 6, dpi = 300)
cat("[保存] 11_stemness_by_cluster.pdf/png\n")

# --- 6. 保存 ---
saveRDS(srt, "03_results/09_single_cell/12_srt_with_stemness.rds")
cat("\n[保存] 12_srt_with_stemness.rds\n")

cat("\n🎉 StemID修复版完成！\n")
