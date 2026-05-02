#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     简化版线粒体模块分析 (Seurat原生 + 相关性网络)         ║\n")
cat("║     替代hdWGCNA，不依赖GitHub安装                        ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(reshape2)
})

srt <- readRDS("03_results/09_single_cell/22_srt_with_pseudotime.rds")
cat("✅ 加载对象:", ncol(srt), "cells\n")

# ============================================================
# 1. 定义线粒体基因集（同hdWGCNA方案）
# ============================================================
cat("\n========== 1. 定义线粒体基因集 ==========\n")

cI_n <- c("NDUFS1","NDUFS2","NDUFS3","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NDUFAB1")
cI_q <- c("NDUFA1","NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10",
          "NDUFA11","NDUFA12","NDUFA13","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5",
          "NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11")
cI_core <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6")
cII <- c("SDHA","SDHB","SDHC","SDHD")
cIII <- c("UQCRC1","UQCRC2","CYC1","UQCRB","UQCRQ","UQCR10","UQCR11")
cIV <- c("COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A1","COX7B","COX8A")
cV <- c("ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5MC1","ATP5MC2","ATP5MC3")
mito_bio <- c("PPARGC1A","NRF1","NRF2","TFAM","TFB1M","TFB2M")
mito_dyn <- c("MFN1","MFN2","OPA1","DNM1L","FIS1")

all_mito <- unique(c(cI_n, cI_q, cI_core, cII, cIII, cIV, cV, mito_bio, mito_dyn))
avail_mito <- all_mito[all_mito %in% rownames(srt)]
cat("线粒体基因集:", length(all_mito), "个定义 |", length(avail_mito), "个可用\n")

# 扩展：加入所有OXPHOS相关基因（更全面的模块搜索）
oxphos_genes <- c(avail_mito, grep("^NDUF", rownames(srt), value = TRUE),
                  grep("^COX", rownames(srt), value = TRUE),
                  grep("^ATP5", rownames(srt), value = TRUE),
                  grep("^UQC", rownames(srt), value = TRUE),
                  grep("^SDH", rownames(srt), value = TRUE))
oxphos_genes <- unique(oxphos_genes)
oxphos_genes <- oxphos_genes[oxphos_genes %in% rownames(srt)]
cat("扩展OXPHOS基因:", length(oxphos_genes), "个\n")

# ============================================================
# 2. NDUFB7共表达邻居分析（模块核心）
# ============================================================
cat("\n========== 2. NDUFB7共表达邻居分析 ==========\n")

expr_mat <- as.matrix(GetAssayData(srt, assay = "RNA", slot = "data")[oxphos_genes, ])

# 计算NDUFB7与所有OXPHOS基因的Spearman相关
ndufb7_expr <- as.numeric(expr_mat["NDUFB7", ])
cor_results <- data.frame(gene = oxphos_genes, rho = NA, p = NA)

for (i in 1:nrow(expr_mat)) {
  g <- oxphos_genes[i]
  if (g != "NDUFB7") {
    ct <- cor.test(ndufb7_expr, as.numeric(expr_mat[g, ]), method = "spearman")
    cor_results$rho[i] <- ct$estimate
    cor_results$p[i] <- ct$p.value
  } else {
    cor_results$rho[i] <- 1.0
    cor_results$p[i] <- 0
  }
}

cor_results <- cor_results[order(-cor_results$rho), ]
cor_results$FDR <- p.adjust(cor_results$p, method = "fdr")
cor_results$module_member <- cor_results$rho > 0.2 & cor_results$FDR < 0.05

cat("Top 20 NDUFB7共表达邻居:\n")
print(head(cor_results[, c("gene","rho","p","FDR","module_member")], 20))

# 模块基因列表
module_genes <- cor_results$gene[cor_results$module_member]
if (length(module_genes) < 5) {
  cat("⚠️ 严格阈值下模块基因不足，放宽至rho>0.15\n")
  module_genes <- cor_results$gene[cor_results$rho > 0.15 & cor_results$p < 0.05]
}
cat("\n🎯 NDUFB7模块基因数:", length(module_genes), "\n")
cat("模块基因:", paste(module_genes, collapse = ", "), "\n")

write.csv(cor_results, "03_results/09_single_cell/23_ndufb7_correlation_network.csv", row.names = FALSE)
cat("[保存] 23_ndufb7_correlation_network.csv\n")

# 检查PGC-1α/NRF1/TFAM是否在模块中
buddies <- c("PPARGC1A","NRF1","NRF2","TFAM")
for (g in buddies) {
  if (g %in% module_genes) {
    cat("✅", g, "在NDUFB7模块中！rho =", round(cor_results$rho[cor_results$gene==g], 3), "\n")
  } else if (g %in% cor_results$gene) {
    cat("⚠️", g, "不在模块中。rho =", round(cor_results$rho[cor_results$gene==g], 3), 
        ", p =", format(cor_results$p[cor_results$gene==g], digits=3), "\n")
  } else {
    cat("❌", g, "不在数据中\n")
  }
}

# ============================================================
# 3. 模块评分计算（AddModuleScore）
# ============================================================
cat("\n========== 3. 模块评分计算 ==========\n")

if (length(module_genes) >= 3) {
  srt <- AddModuleScore(srt, features = list(module_genes), name = "NDUFB7_Module", assay = "RNA")
  srt$ndufb7_module_score <- srt$NDUFB7_Module1
  
  cat("模块评分范围:", round(range(srt$ndufb7_module_score, na.rm=TRUE), 3), "\n")
  
  # 按Cluster统计
  mod_cluster <- srt@meta.data %>%
    group_by(seurat_clusters, condition) %>%
    summarise(
      mean_score = mean(ndufb7_module_score, na.rm = TRUE),
      median_score = median(ndufb7_module_score, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  cat("\n模块评分按Cluster:\n")
  print(as.data.frame(mod_cluster))
  
  # Cluster 4 vs 其他HF的检验
  c4_score <- srt$ndufb7_module_score[srt$seurat_clusters == 4]
  other_hf_score <- srt$ndufb7_module_score[srt$condition == "HF" & srt$seurat_clusters != 4]
  wt <- wilcox.test(c4_score, other_hf_score)
  cat("\nCluster 4 vs 其他HF模块评分: p =", format(wt$p.value, digits=3), "\n")
  cat("Cluster 4均值:", round(mean(c4_score), 3), "| 其他HF均值:", round(mean(other_hf_score), 3), "\n")
  
  # HF vs Control
  hf_score <- srt$ndufb7_module_score[srt$condition == "HF"]
  ctrl_score <- srt$ndufb7_module_score[srt$condition == "Control"]
  wt2 <- wilcox.test(hf_score, ctrl_score)
  cat("HF vs Control模块评分: p =", format(wt2$p.value, digits=3), "\n")
  cat("HF均值:", round(mean(hf_score), 3), "| Control均值:", round(mean(ctrl_score), 3), "\n")
  
} else {
  cat("❌ 模块基因不足，跳过模块评分\n")
}

# ============================================================
# 4. 可视化
# ============================================================
cat("\n========== 4. 可视化 ==========\n")

# 4A: 相关性网络热图（Top 30邻居）
top_n <- min(30, nrow(cor_results))
top_genes <- head(cor_results$gene, top_n)
cor_sub <- cor(t(expr_mat[top_genes, ]), method = "spearman")

pdf("03_results/09_single_cell/30_ndufb7_module_heatmap.pdf", width = 10, height = 9)
pheatmap(cor_sub, 
         display_numbers = TRUE, 
         number_format = "%.2f",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = paste("NDUFB7 Co-expression Module (Top", top_n, "genes)"),
         fontsize_number = 6)
dev.off()
cat("[保存] 30_ndufb7_module_heatmap.pdf\n")

# 4B: 模块评分UMAP
p_mod <- FeaturePlot(srt, features = "ndufb7_module_score", pt.size = 1.5, order = TRUE) +
  scale_color_gradientn(colors = c("navy", "cyan", "yellow", "red"), name = "Module Score") +
  ggtitle("NDUFB7 Module Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("03_results/09_single_cell/31_module_score_umap.pdf", p_mod, width = 8, height = 6)
ggsave("03_results/09_single_cell/31_module_score_umap.png", p_mod, width = 8, height = 6, dpi = 300)
cat("[保存] 31_module_score_umap.pdf/png\n")

# 4C: 模块评分Violin（按Cluster+Condition）
p_vln <- VlnPlot(srt, features = "ndufb7_module_score", pt.size = 0.3, split.by = "condition") +
  ggtitle("NDUFB7 Module Score by Cluster") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("03_results/09_single_cell/32_module_score_violin.pdf", p_vln, width = 12, height = 6)
ggsave("03_results/09_single_cell/32_module_score_violin.png", p_vln, width = 12, height = 6, dpi = 300)
cat("[保存] 32_module_score_violin.pdf/png\n")

# 4D: 模块基因在Cluster中的表达热图
if (length(module_genes) >= 3) {
  avg_mod <- AverageExpression(srt, features = module_genes, group.by = "seurat_clusters")$RNA
  pdf("03_results/09_single_cell/33_module_genes_cluster_heatmap.pdf", width = 8, height = max(6, length(module_genes)*0.3))
  pheatmap(avg_mod, scale = "row",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "NDUFB7 Module Genes by Cluster (Z-score)")
  dev.off()
  cat("[保存] 33_module_genes_cluster_heatmap.pdf\n")
}

# ============================================================
# 5. 保存
# ============================================================
cat("\n========== 5. 保存 ==========\n")
saveRDS(srt, "03_results/09_single_cell/34_srt_with_module.rds")
cat("[保存] 34_srt_with_module.rds\n")

cat("\n🎉 简化版模块分析完成！\n")
cat("\n【关键产出】\n")
cat("  23_ndufb7_correlation_network.csv —— NDUFB7与所有OXPHOS基因的相关性\n")
cat("  30_ndufb7_module_heatmap.pdf —— 模块相关性热图\n")
cat("  31_module_score_umap.pdf —— 模块评分UMAP\n")
cat("  32_module_score_violin.pdf —— 模块评分分布\n")
cat("  33_module_genes_cluster_heatmap.pdf —— 模块基因Cluster热图\n")
