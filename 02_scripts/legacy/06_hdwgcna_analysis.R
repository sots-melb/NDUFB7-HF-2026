#!/usr/bin/env Rscript
# 02_scripts/05_scRNA_analysis/06_hdwgcna_analysis.R

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     hdWGCNA 线粒体共表达模块分析 (Seurat V4 环境)          ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

library(Seurat)
library(hdWGCNA)
library(dplyr)
library(ggplot2)
library(patchwork)
library(WGCNA)

# 启用WGCNA多线程以加速计算
allowWGCNAThreads()

PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)
dir.create("03_results/09_single_cell/hdwgcna", recursive = TRUE, showWarnings = FALSE)

cat("========== 1. 读取单细胞数据 ==========\n")
# 读取包含Cluster信息和简化模块的Seurat对象
srt <- readRDS("03_results/09_single_cell/34_srt_with_module.rds")

cat("========== 2. hdWGCNA 初始化 ==========\n")
# SetupForWGCNA: 选择基因进行分析
srt <- SetupForWGCNA(
  srt,
  gene_select = "fraction", # 选择在一定比例细胞中表达的基因
  fraction = 0.05,          # 至少在5%的细胞中表达
  wgcna_name = "CM_OXPHOS"
)

cat("========== 3. 构建 Metacells (元细胞) ==========\n")
# 克服scRNA-seq的稀疏性
# 【终极修复】：移除 condition 分组，避免极端稀疏的子集（如某组仅1个细胞）导致 KNN 计算崩溃
srt <- MetacellsByGroups(
  seurat_obj = srt,
  group.by = "seurat_clusters", # 仅按 Cluster 聚合，不再按 HF/Control 切碎
  k = 15,                # KNN近邻数
  max_shared = 5,        # 允许最大共享细胞数
  min_cells = 20,        # 最小细胞数要求：20，足以保住 32 个细胞的 Cluster 4
  ident.group = 'seurat_clusters'
)
srt <- NormalizeMetacells(srt)

cat("========== 4. 软阈值 (Soft Power) 测试 ==========\n")
# 修复保存报错：使用基础的 pdf() 打印列表对象
pdf("03_results/09_single_cell/hdwgcna/01_soft_power.pdf", width = 10, height = 8)
print(p_power)
dev.off()
cat("✅ 软阈值图已保存\n")

cat("========== 5. 构建共表达网络 ==========\n")
# ConstructNetwork 会自动使用上面计算出的最佳软阈值
srt <- ConstructNetwork(
  srt,
  tom_name = 'CM_OXPHOS',
  overwrite_tom = TRUE
)

# 绘制基因树状图
pdf("03_results/09_single_cell/hdwgcna/02_dendrogram.pdf", width=8, height=6)
PlotDendrogram(srt, main='hdWGCNA Dendrogram (CM OXPHOS)')
dev.off()
cat("✅ 基因树状图已保存\n")

cat("========== 6. 计算模块特征基因 (Module Eigengenes) ==========\n")
srt <- ModuleEigengenes(srt)
srt <- ModuleConnectivity(srt)

cat("========== 7. 定位 NDUFB7 所在模块 ==========\n")
modules <- GetModules(srt)
ndufb7_module <- modules %>% filter(gene_name == "NDUFB7") %>% pull(module)

if(length(ndufb7_module) > 0) {
  cat("🎯 核心发现: NDUFB7 被分配在模块 [", ndufb7_module, "]\n")
  
  # 1. 保存 Hub Gene 网络图
  pdf(paste0("03_results/09_single_cell/hdwgcna/03_hub_network_", ndufb7_module, ".pdf"), width = 8, height = 8)
  p_hub <- HubGeneNetworkPlot(srt, target_modules = ndufb7_module, n_hubs = 10, n_other = 5, edge_alpha = 0.5)
  print(p_hub)
  dev.off()
  
  # 2. 模块特征分布图 (验证 Cluster 4 塌陷)
  module_name <- paste0("CM_OXPHOS-", ndufb7_module)
  if(module_name %in% colnames(srt@meta.data)) {
    pdf(paste0("03_results/09_single_cell/hdwgcna/04_module_vln_", ndufb7_module, ".pdf"), width = 10, height = 6)
    p_vln <- VlnPlot(srt, features = module_name, group.by = "seurat_clusters", split.by = "condition", pt.size = 0) +
      ggtitle(paste("Module Eigengene:", ndufb7_module, "(Contains NDUFB7)"))
    print(p_vln)
    dev.off()
  }
} else {
  cat("⚠️ 警告: NDUFB7 未被分配到任何模块\n")
}

cat("========== 8. 保存最终 hdWGCNA 对象 ==========\n")
saveRDS(srt, "03_results/09_single_cell/hdwgcna/srt_hdwgcna_completed.rds")
cat("🎉 hdWGCNA 完整分析顺利结束！\n")