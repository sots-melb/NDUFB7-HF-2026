#!/usr/bin/env Rscript
# GSE315590 小鼠TAC跨物种验证 - 完整版
setwd("~/Projects/NDUFB7_HF_2026_04_20")

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     GSE315590 小鼠TAC跨物种验证                          ║\n")
cat("║     目标：Ndufb7在压力负荷心衰中的表达模式               ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# --- 0. 包加载 ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

out_dir <- "03_results/20_mouse_tac"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. 读取数据 ---
counts_file <- "01_data/01_raw_geo/GSE315590/scRNAseq_TAC_transition_raw_counts.txt.gz"
if (!file.exists(counts_file)) {
  cat("❌ 数据文件不存在:", counts_file, "\n")
  quit(status = 1)
}

fsize <- round(file.info(counts_file)$size / 1024 / 1024, 2)
cat("✅ 数据文件存在:", fsize, "MB\n")

cat("\n========== 1. 读取小鼠TAC数据 ==========\n")
counts <- read.table(gzfile(counts_file), header = TRUE, row.names = 1, 
                     sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
cat("原始矩阵:", nrow(counts), "genes x", ncol(counts), "cells\n")

# --- 2. 检查Ndufb7 ---
mouse_gene <- "Ndufb7"
if (mouse_gene %in% rownames(counts)) {
  cat("✅ 找到小鼠同源基因:", mouse_gene, "\n")
} else if ("Ndufb7" %in% rownames(counts)) {
  mouse_gene <- "Ndufb7"
  cat("✅ 找到小鼠同源基因(首字母大写):", mouse_gene, "\n")
} else {
  cat("⚠️  未找到Ndufb7，尝试模糊匹配...\n")
  candidates <- grep("ndufb7", rownames(counts), ignore.case = TRUE, value = TRUE)
  if (length(candidates) > 0) {
    mouse_gene <- candidates[1]
    cat("✅ 模糊匹配到:", mouse_gene, "\n")
  } else {
    cat("❌ 完全未找到Ndufb7同源基因\n")
    quit(status = 1)
  }
}

# --- 3. 创建Seurat对象 ---
cat("\n========== 2. 创建Seurat对象 ==========\n")
# 假设列名格式包含分组信息 (如 Sham_1, TAC_1 等)
# 尝试从列名推断分组
meta <- data.frame(cell = colnames(counts), row.names = colnames(counts))

# 推断分组：如果列名包含Sham/TAC/Sham_Week/TAC_Week等
if (any(grepl("Sham|TAC", colnames(counts), ignore.case = TRUE))) {
  meta$group <- ifelse(grepl("Sham", colnames(counts), ignore.case = TRUE), "Sham", "TAC")
  # 尝试推断时间点
  if (any(grepl("week|day|w\\d|d\\d", colnames(counts), ignore.case = TRUE))) {
    time_match <- regmatches(colnames(counts), 
      regexpr("(?i)(week[_]?\\d+|day[_]?\\d+|w\\d+|d\\d+)", colnames(counts)))
    meta$time <- ifelse(length(time_match) == ncol(counts), time_match, "unknown")
  } else {
    meta$time <- "unknown"
  }
} else {
  meta$group <- "unknown"
  meta$time <- "unknown"
}

cat("分组推断:\n")
print(table(meta$group))

seu <- CreateSeuratObject(counts = counts, meta.data = meta, project = "Mouse_TAC")
cat("Seurat对象创建完成:", ncol(seu), "cells\n")

# --- 4. QC ---
cat("\n========== 3. 质控过滤 ==========\n")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & 
              nCount_RNA > 500 & percent.mt < 15)
cat("QC后:", ncol(seu), "cells\n")

# --- 5. 标准化 ---
cat("\n========== 4. 标准化与降维 ==========\n")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.8)
seu <- RunUMAP(seu, dims = 1:20)

cat("聚类完成，共", length(unique(seu$seurat_clusters)), "个cluster\n")

# --- 6. Ndufb7表达可视化 ---
cat("\n========== 5. Ndufb7表达分析 ==========\n")

# UMAP图
p1 <- DimPlot(seu, reduction = "umap", label = TRUE) + ggtitle("Mouse TAC Clusters")
p2 <- FeaturePlot(seu, features = mouse_gene, order = TRUE, 
                  cols = c("grey", "red")) + ggtitle(paste(mouse_gene, "Expression"))

ggsave(file.path(out_dir, "01_umap_clusters.png"), p1, width = 8, height = 6, dpi = 150)
ggsave(file.path(out_dir, "02_ndufb7_umap.png"), p2, width = 8, height = 6, dpi = 150)
cat("✅ UMAP图已保存\n")

# 小提琴图
p3 <- VlnPlot(seu, features = mouse_gene, pt.size = 0) + 
  ggtitle(paste(mouse_gene, "by Cluster"))
ggsave(file.path(out_dir, "03_ndufb7_vln.png"), p3, width = 10, height = 6, dpi = 150)

# 按分组比较
if ("group" %in% colnames(seu@meta.data) && length(unique(seu$group)) > 1) {
  p4 <- VlnPlot(seu, features = mouse_gene, group.by = "group", pt.size = 0) +
    ggtitle(paste(mouse_gene, "by Group (Sham vs TAC)"))
  ggsave(file.path(out_dir, "04_ndufb7_by_group.png"), p4, width = 8, height = 6, dpi = 150)
  
  # 统计检验
  sham_expr <- GetAssayData(seu, layer = "data")[mouse_gene, seu$group == "Sham"]
  tac_expr <- GetAssayData(seu, layer = "data")[mouse_gene, seu$group == "TAC"]
  wilcox_res <- wilcox.test(sham_expr, tac_expr)
  cat("Sham vs TAC Wilcoxon p-value:", format(wilcox_res$p.value, digits = 3), "\n")
}

# --- 7. 各Cluster Ndufb7表达量统计 ---
cat("\n========== 6. Cluster-level统计 ==========\n")
cluster_expr <- AverageExpression(seu, features = mouse_gene, group.by = "seurat_clusters")
cluster_expr_df <- as.data.frame(cluster_expr$RNA)
cluster_expr_df$cluster <- rownames(cluster_expr_df)
colnames(cluster_expr_df)[1] <- "avg_expression"
cluster_expr_df <- cluster_expr_df[order(-cluster_expr_df$avg_expression), ]
write.csv(cluster_expr_df, file.path(out_dir, "05_ndufb7_cluster_avg.csv"), row.names = FALSE)
cat("Cluster平均表达:\n")
print(head(cluster_expr_df, 10))

# 找出高表达cluster
high_cluster <- cluster_expr_df$cluster[1]
cat("\n🏆 Ndufb7最高表达Cluster:", high_cluster, 
    "(avg=", round(cluster_expr_df$avg_expression[1], 4), ")\n")

# --- 8. 跨物种对比注释 ---
cat("\n========== 7. 跨物种注释推断 ==========\n")
# 基于marker基因推断cluster身份
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, file.path(out_dir, "06_top10_markers.csv"), row.names = FALSE)

# 简单注释规则
cluster_anno <- data.frame(
  cluster = as.character(0:(length(unique(seu$seurat_clusters))-1)),
  annotation = "Unknown"
)
# 基于常见心脏细胞marker进行推断
for (cl in unique(markers$cluster)) {
  cl_markers <- markers$gene[markers$cluster == cl]
  if (any(grepl("Myh6|Myh7|Actc1|Ttn", cl_markers, ignore.case = TRUE))) {
    cluster_anno$annotation[cluster_anno$cluster == cl] <- "Cardiomyocyte"
  } else if (any(grepl("Pecam1|Cdh5|Vwf", cl_markers, ignore.case = TRUE))) {
    cluster_anno$annotation[cluster_anno$cluster == cl] <- "Endothelial"
  } else if (any(grepl("Col1a1|Col3a1|Postn", cl_markers, ignore.case = TRUE))) {
    cluster_anno$annotation[cluster_anno$cluster == cl] <- "Fibroblast"
  } else if (any(grepl("Ptprc|Cd68|Cd14", cl_markers, ignore.case = TRUE))) {
    cluster_anno$annotation[cluster_anno$cluster == cl] <- "Immune"
  }
}
write.csv(cluster_anno, file.path(out_dir, "07_cluster_annotation.csv"), row.names = FALSE)
cat("Cluster注释已保存\n")

# --- 9. 保存Seurat对象 ---
cat("\n========== 8. 保存结果 ==========\n")
saveRDS(seu, file.path(out_dir, "mouse_tac_seurat.rds"))
cat("✅ Seurat对象已保存:", file.path(out_dir, "mouse_tac_seurat.rds"), "\n")

# 保存Ndufb7表达矩阵
ndufb7_expr <- data.frame(
  cell = colnames(seu),
  group = seu$group,
  cluster = seu$seurat_clusters,
  ndufb7 = GetAssayData(seu, layer = "data")[mouse_gene, ]
)
write.csv(ndufb7_expr, file.path(out_dir, "08_ndufb7_expression.csv"), row.names = FALSE)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║  小鼠TAC分析完成！                                         ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")
cat("\n📁 结果目录:", out_dir, "\n")
cat("📊 关键结果:\n")
cat("   01_umap_clusters.png       - UMAP聚类图\n")
cat("   02_ndufb7_umap.png         - Ndufb7表达UMAP\n")
cat("   03_ndufb7_vln.png          - Ndufb7小提琴图\n")
cat("   04_ndufb7_by_group.png     - Sham vs TAC比较\n")
cat("   05_ndufb7_cluster_avg.csv  - 各Cluster平均表达\n")
cat("   06_top10_markers.csv       - Top10 marker基因\n")
cat("   07_cluster_annotation.csv  - Cluster注释\n")
cat("   08_ndufb7_expression.csv   - 单细胞表达矩阵\n")
cat("   mouse_tac_seurat.rds       - Seurat对象\n")
