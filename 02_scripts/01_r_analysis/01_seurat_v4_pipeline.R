#!/usr/bin/env Rscript
# NDUFB7心衰研究 - Phase 2 Pillar 2
# GSE168742 单细胞Seurat V4分析
# 环境: 系统默认R 4.3.1 (Seurat 4.3.0.1)

PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     Phase 2 Pillar 2: GSE168742 Seurat V4分析            ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# ==================== 0. 环境检查 ====================
cat("\n========== 0. 环境检查 ==========\n")

library(Seurat)
cat("✅ Seurat版本:", as.character(packageVersion("Seurat")), "\n")

required_pkgs <- c("dplyr", "ggplot2", "patchwork")
for (pkg in required_pkgs) {
  library(pkg, character.only = TRUE)
  cat("✅", pkg, "\n")
}

dir.create("03_results/09_single_cell", showWarnings = FALSE, recursive = TRUE)

# ==================== 1. 读取GSE168742数据 ====================
cat("\n========== 1. 读取GSE168742数据 ==========\n")

# 注意：GSE168742有两个文件：HF和Control
hf_file <- "01_data/01_raw_geo/GSE168742/GSE168742_human_HF_CM_matrix.rds"
ctrl_file <- "01_data/01_raw_geo/GSE168742/GSE168742_human_control_CM.rds"

if (!file.exists(hf_file)) stop("❌ HF文件不存在")
if (!file.exists(ctrl_file)) stop("❌ Control文件不存在")

# 读取HF
obj_hf <- readRDS(hf_file)
cat("HF RDS结构:\n")
print(str(obj_hf))

# 读取Control
obj_ctrl <- readRDS(ctrl_file)
cat("\nControl RDS结构:\n")
print(str(obj_ctrl))

# 提取表达矩阵（适配多种RDS结构）
extract_matrix <- function(obj, label) {
  if (is.matrix(obj) || inherits(obj, "dgCMatrix") || inherits(obj, "dgTMatrix")) {
    mat <- as(obj, "sparseMatrix")
    cat("[", label, "] 读取为稀疏矩阵:", nrow(mat), "x", ncol(mat), "\n")
    return(mat)
  } else if (is.list(obj) && "exprs" %in% names(obj)) {
    mat <- as(obj$exprs, "sparseMatrix")
    cat("[", label, "] 从list提取exprs:", nrow(mat), "x", ncol(mat), "\n")
    return(mat)
  } else {
    stop(paste("❌", label, "未知RDS结构"))
  }
}

exprs_hf <- extract_matrix(obj_hf, "HF")
exprs_ctrl <- extract_matrix(obj_ctrl, "Control")

cat("\nHF矩阵: ", nrow(exprs_hf), "genes x", ncol(exprs_hf), "cells\n")
cat("Control矩阵: ", nrow(exprs_ctrl), "genes x", ncol(exprs_ctrl), "cells\n")

# 检查NDUFB7
cat("\nHF NDUFB7:", "NDUFB7" %in% rownames(exprs_hf), "\n")
cat("Control NDUFB7:", "NDUFB7" %in% rownames(exprs_ctrl), "\n")

# 探测数据类型（counts vs normalized）
sample_vals <- as.numeric(exprs_hf[1:100, 1:10])
cat("\nHF数值示例:\n")
print(head(sample_vals, 20))
cat("整数比例:", round(sum(sample_vals == floor(sample_vals)) / length(sample_vals) * 100, 1), "%\n")
cat("最大值:", max(sample_vals), "\n")

# 如果已经是log-normalized（非整数、值较小），需要反向转换或标记
is_log_normalized <- max(sample_vals) < 100 && mean(sample_vals) < 10
if (is_log_normalized) {
  cat("⚠️ 数据疑似已log-normalized，将标记为data slot而非counts\n")
}

# ==================== 2. 创建Seurat对象 ====================
cat("\n========== 2. 创建Seurat V4对象 ==========\n")

# 分别创建HF和Control对象
srt_hf <- CreateSeuratObject(counts = exprs_hf, project = "HF", min.cells = 3, min.features = 200)
srt_ctrl <- CreateSeuratObject(counts = exprs_ctrl, project = "Control", min.cells = 3, min.features = 200)

# 添加condition标签
srt_hf$condition <- "HF"
srt_ctrl$condition <- "Control"

cat("HF对象:", ncol(srt_hf), "cells,", nrow(srt_hf), "genes\n")
cat("Control对象:", ncol(srt_ctrl), "cells,", nrow(srt_ctrl), "genes\n")

# 合并（V4语法）
srt <- merge(srt_hf, srt_ctrl)
cat("\n合并后:", ncol(srt), "cells,", nrow(srt), "genes\n")
cat("Condition分布:\n")
print(table(srt$condition))

# 如果数据是log-normalized，需要特殊处理
if (is_log_normalized) {
  cat("\n⚠️ 使用预计算数据模式（非原始counts）\n")
  # 将数据存入data slot，跳过NormalizeData
  srt@assays$RNA@data <- srt@assays$RNA@counts
}

# ==================== 3. QC过滤 ====================
cat("\n========== 3. QC过滤 ==========\n")

# 计算线粒体比例
srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^MT-")

cat("QC前统计:\n")
cat("  总细胞:", ncol(srt), "\n")
cat("  nFeature_RNA:", paste(round(quantile(srt$nFeature_RNA), 0), collapse = " / "), "\n")
cat("  nCount_RNA:", paste(round(quantile(srt$nCount_RNA), 0), collapse = " / "), "\n")
cat("  percent.mt:", paste(round(quantile(srt$percent.mt), 2), collapse = " / "), "%\n")

# 保存QC图
pdf("03_results/09_single_cell/01_QC_violin.pdf", width = 12, height = 5)
VlnPlot(srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0.1, group.by = "condition")
dev.off()
cat("[保存] 01_QC_violin.pdf\n")

# 过滤（心肌细胞参数）
srt <- subset(srt, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)

cat("\nQC后:\n")
cat("  剩余细胞:", ncol(srt), "\n")
cat("  过滤比例:", round((1 - ncol(srt)/(ncol(srt_hf)+ncol(srt_ctrl)))*100, 1), "%\n")
cat("  Condition分布:\n")
print(table(srt$condition))

# ==================== 4. 标准化与特征选择 ====================
cat("\n========== 4. 标准化与特征选择 ==========\n")

if (!is_log_normalized) {
  # 原始counts，正常标准化
  srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
  cat("✅ LogNormalize完成\n")
} else {
  cat("⏭️ 跳过NormalizeData（数据已标准化）\n")
}

srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)
cat("高变基因:", length(VariableFeatures(srt)), "\n")

# 缩放
srt <- ScaleData(srt)
cat("✅ ScaleData完成\n")

# ==================== 5. 降维与聚类 ====================
cat("\n========== 5. 降维与聚类 ==========\n")

srt <- RunPCA(srt, features = VariableFeatures(object = srt))
cat("PCA完成, 维度:", length(srt@reductions$pca@cell.embeddings[1,]), "\n")

# 选择PC（使用ElbowPlot或自动选择）
pdf("03_results/09_single_cell/02_PCA_elbow.pdf", width = 8, height = 5)
ElbowPlot(srt, ndims = 30)
dev.off()
cat("[保存] 02_PCA_elbow.pdf\n")

# 使用20个PC
ndims <- 20
srt <- FindNeighbors(srt, dims = 1:ndims)
srt <- FindClusters(srt, resolution = 0.8)
srt <- RunUMAP(srt, dims = 1:ndims)

cat("聚类数:", length(levels(srt$seurat_clusters)), "\n")
cat("每簇细胞数:\n")
print(table(srt$seurat_clusters, srt$condition))

# ==================== 6. 可视化 ====================
cat("\n========== 6. 可视化 ==========\n")

# UMAP分condition
pdf("03_results/09_single_cell/03_UMAP_condition.pdf", width = 10, height = 5)
p1 <- DimPlot(srt, reduction = "umap", group.by = "condition", pt.size = 0.5) + 
      ggtitle("Condition")
p2 <- DimPlot(srt, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) + 
      ggtitle("Clusters")
print(p1 + p2)
dev.off()
cat("[保存] 03_UMAP_condition.pdf\n")

# ==================== 7. NDUFB7表达分析 ====================
cat("\n========== 7. NDUFB7表达分析 ==========\n")

if ("NDUFB7" %in% rownames(srt)) {
  cat("✅ NDUFB7在数据中\n")
  
  # FeaturePlot + VlnPlot
  pdf("03_results/09_single_cell/04_NDUFB7_expression.pdf", width = 12, height = 5)
  p1 <- FeaturePlot(srt, features = "NDUFB7", pt.size = 0.5, cols = c("grey", "red")) +
        ggtitle("NDUFB7 Expression")
  p2 <- VlnPlot(srt, features = "NDUFB7", pt.size = 0.1, group.by = "condition") +
        ggtitle("NDUFB7 by Condition")
  print(p1 + p2)
  dev.off()
  cat("[保存] 04_NDUFB7_expression.pdf\n")
  
  # 各簇NDUFB7统计
  ndufb7_data <- FetchData(srt, vars = c("NDUFB7", "seurat_clusters", "condition"))
  cluster_stats <- aggregate(NDUFB7 ~ seurat_clusters + condition, data = ndufb7_data, FUN = mean)
  cat("\n各簇-Condition NDUFB7平均表达:\n")
  print(cluster_stats[order(-cluster_stats$NDUFB7), ])
  
  # 保存统计
  write.csv(cluster_stats, "03_results/09_single_cell/05_NDUFB7_cluster_stats.csv", row.names = FALSE)
  cat("[保存] 05_NDUFB7_cluster_stats.csv\n")
  
} else {
  cat("❌ NDUFB7不在数据中，检查行名...\n")
  hits <- grep("NDUFB7|67996", rownames(srt), value = TRUE, ignore.case = TRUE)
  cat("匹配:", paste(head(hits, 5), collapse = ", "), "\n")
}

# ==================== 8. 保存对象 ====================
cat("\n========== 8. 保存Seurat对象 ==========\n")

saveRDS(srt, "03_results/09_single_cell/06_srt_v4_processed.rds")
cat("[保存] 06_srt_v4_processed.rds (", round(file.info("03_results/09_single_cell/06_srt_v4_processed.rds")$size/1024/1024, 1), "MB)\n")

cat("\n========================================\n")
cat("🎉 Seurat V4基础分析完成！\n")
cat("========================================\n")
cat("\n下一步:\n")
cat("1. 检查UMAP和NDUFB7表达图\n")
cat("2. 手动判断哪些cluster是Cardiomyocyte\n")
cat("3. 提取CM子集进行StemID分析\n")
