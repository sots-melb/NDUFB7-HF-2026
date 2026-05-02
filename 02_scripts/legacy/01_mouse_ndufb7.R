#!/usr/bin/env Rscript
setwd("~/Projects/NDUFB7_HF_{2026_04_20}")

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     GSE315590 小鼠TAC跨物种验证                          ║\n")
cat("║     目标：验证Cluster 4型亚群是否在小鼠压力负荷心衰复现  ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# --- 0. 检查数据 ---
counts_file <- "01_data/01_raw_geo/GSE315590/scRNAseq_TAC_transition_raw_counts.txt.gz"
if (!file.exists(counts_file)) {
  cat("❌ 数据文件不存在:", counts_file, "\n")
  cat("请确认GSE315590数据已下载到 01_data/01_raw_geo/GSE315590/\n")
  quit(status = 1)
}

cat("✅ 数据文件存在，大小:", round(file.info(counts_file)$size/1024/1024, 2), "MB\n")

# --- 1. 读取数据 ---
cat("\n========== 1. 读取小鼠TAC数据 ==========\n")
# 读取count矩阵（可能是gene x cell格式）
counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("原始矩阵:", nrow(counts), "genes x", ncol(counts), "cells\n")

# 检查Ndufb7（小鼠同源基因，首字母小写）
mouse_ndufb7 <- "Ndufb7"
if (mouse_ndufb7 %in% rownames(counts)) {
  cat("✅ 找到小鼠同源基因:", mouse_ndufb7, "\n")
} else if ("Ndufb7" %in% rownames(counts)) {
  cat("✅ 找到小鼠同源基因（大小写匹配）:", "Ndufb7", "\n")
  mouse_ndufb7 <- "Ndufb7"
} else {
  cat("⚠️ 未找到Ndufb7，尝试搜索NDUFB7...\n")
  candidates <- grep("ndufb7", rownames(counts), ignore.case = TRUE, value = TRUE)
  cat("候选:", paste(candidates, collapse = ", "), "\n")
  if (length(candidates) > 0) mouse_ndufb7 <- candidates[1]
}

# --- 2. 构建Seurat对象 ---
cat("\n========== 2. 构建Seurat对象 ==========\n")
srt_mouse <- CreateSeuratObject(counts = counts, project = "TAC", min.cells = 3, min.features = 200)
cat("Seurat对象:", ncol(srt_mouse), "cells x", nrow(srt_mouse), "genes\n")

# 从列名推断时间点和样本
# GSE315590通常是TAC时间序列（sham, 3d, 7d, 14d, 28d等）
colnames_orig <- colnames(srt_mouse)
cat("前10个细胞名:", paste(head(colnames_orig, 10), collapse = ", "), "\n")

# 尝试从细胞名提取时间信息（常见格式: sample_time_barcode）
# 如果无法自动推断，需要手动注释
time_pattern <- grep("sham|day|d[0-9]|week|w[0-9]", colnames_orig, ignore.case = TRUE, value = TRUE)
cat("检测到时间相关命名:", length(time_pattern), "个\n")

# 简化：假设用户需要手动提供分组信息，或从GEO metadata获取
# 这里先进行基础QC和聚类
srt_mouse[["percent.mt"]] <- PercentageFeatureSet(srt_mouse, pattern = "^mt-")
srt_mouse <- subset(srt_mouse, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 15)
cat("QC后:", ncol(srt_mouse), "cells\n")

srt_mouse <- NormalizeData(srt_mouse)
srt_mouse <- FindVariableFeatures(srt_mouse, selection.method = "vst", nfeatures = 2000)
srt_mouse <- ScaleData(srt_mouse)
srt_mouse <- RunPCA(srt_mouse, features = VariableFeatures(object = srt_mouse))
srt_mouse <- FindNeighbors(srt_mouse, dims = 1:15)
srt_mouse <- FindClusters(srt_mouse, resolution = 0.6)
srt_mouse <- RunUMAP(srt_mouse, dims = 1:15)

cat("✅ 聚类完成:", length(unique(srt_mouse$seurat_clusters)), "clusters\n")

# --- 3. Ndufb7分析 ---
cat("\n========== 3. Ndufb7表达分析 ==========\n")
if (mouse_ndufb7 %in% rownames(srt_mouse)) {
  # FeaturePlot
  p1 <- FeaturePlot(srt_mouse, features = mouse_ndufb7, cols = c("grey", "red"), pt.size = 0.5) +
    ggtitle(paste("Ndufb7 (mouse homolog) in TAC"))
  ggsave("03_results/20_mouse_tac/60_mouse_ndufb7_umap.pdf", p1, width = 6, height = 5)
  cat("[保存] 60_mouse_ndufb7_umap.pdf\n")
  
  # Violin by cluster
  p2 <- VlnPlot(srt_mouse, features = mouse_ndufb7, pt.size = 0) +
    ggtitle("Ndufb7 by Cluster (TAC)")
  ggsave("03_results/20_mouse_tac/61_mouse_ndufb7_violin.pdf", p2, width = 6, height = 4)
  cat("[保存] 61_mouse_ndufb7_violin.pdf\n")
  
  # 统计各cluster的Ndufb7
  ndufb7_stats <- FetchData(srt_mouse, vars = c(mouse_ndufb7, "seurat_clusters")) %>%
    group_by(seurat_clusters) %>%
    summarise(mean = mean(!!sym(mouse_ndufb7)), median = median(!!sym(mouse_ndufb7)), n = n()) %>%
    arrange(mean)
  cat("\nNdufb7各Cluster统计:\n")
  print(ndufb7_stats, n = Inf)
  
  # 检查是否有极低表达cluster（类似人Cluster 4）
  min_mean <- min(ndufb7_stats$mean)
  if (min_mean < 0.3) {
    cat("\n🎯 发现小鼠极低表达Cluster（均值<0.3），与人Cluster 4特征一致！\n")
  }
  
  # 保存
  saveRDS(srt_mouse, "03_results/20_mouse_tac/62_srt_mouse_tac.rds")
  write.csv(ndufb7_stats, "03_results/20_mouse_tac/63_mouse_ndufb7_stats.csv", row.names = FALSE)
  cat("\n✅ 小鼠TAC分析完成\n")
  
} else {
  cat("❌ 无法找到Ndufb7同源基因，跳过分析\n")
}

cat("\n【说明】小鼠TAC数据需要手动添加condition/time注释才能进行HF vs Control比较。\n")
cat("建议: 从GEO下载metadata，或根据细胞命名规则推断分组。\n")
