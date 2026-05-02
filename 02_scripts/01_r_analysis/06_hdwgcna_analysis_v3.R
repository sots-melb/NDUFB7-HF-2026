#!/usr/bin/env Rscript
library(Seurat)
library(hdWGCNA)
library(dplyr)
library(ggplot2)
library(patchwork)
library(WGCNA)

# 【核心修复1】：严格限制线程数为 16，防止 CPU 调度雪崩
allowWGCNAThreads(nThreads = 16)
setwd("~/Projects/NDUFB7_HF_{2026_04_20}")

cat("\n========== 1. 读取数据 ==========\n")
srt <- readRDS("03_results/09_single_cell/34_srt_with_module.rds")

cat("\n========== 2. 初始化 ==========\n")
srt <- SetupForWGCNA(srt, gene_select = "fraction", fraction = 0.05, wgcna_name = "CM_OXPHOS") 

cat("\n========== 3. 构建 Metacells ==========\n")
srt <- MetacellsByGroups(
  seurat_obj = srt,
  group.by = "seurat_clusters", 
  k = 15, max_shared = 5, min_cells = 20, 
  ident.group = 'seurat_clusters'
)
srt <- NormalizeMetacells(srt)

cat("\n========== 4. 软阈值测试 ==========\n")
metacell_obj <- GetMetacellObject(srt)
valid_groups <- as.character(unique(metacell_obj$seurat_clusters))
srt <- SetDatExpr(srt, group_name = valid_groups, group.by = 'seurat_clusters')

srt <- TestSoftPowers(srt, networkType = 'signed')
p_power <- PlotSoftPowers(srt)

pdf("03_results/09_single_cell/hdwgcna/01_soft_power.pdf", width = 10, height = 8)
print(p_power)
dev.off()

cat("\n========== 5. 构建网络 ==========\n")
# 【核心修复2】：将 maxBlockSize 设置为 30000 (远大于我们的 15000 个基因)
# 强制 WGCNA 在 128G 内存中一次性算完，绝不写硬盘临时文件，彻底告别 error reading from connection！
srt <- ConstructNetwork(
  srt, 
  tom_name = 'CM_OXPHOS', 
  overwrite_tom = TRUE,
  maxBlockSize = 30000 
)

pdf("03_results/09_single_cell/hdwgcna/02_dendrogram.pdf", width=8, height=6)
PlotDendrogram(srt, main='hdWGCNA Dendrogram')
dev.off()

cat("\n========== 6. 模块特征与导出 ==========\n")
srt <- ModuleEigengenes(srt)
srt <- ModuleConnectivity(srt)

modules <- GetModules(srt)
ndufb7_module <- modules %>% filter(gene_name == "NDUFB7") %>% pull(module)

if(length(ndufb7_module) > 0) {
  cat("\n🎯 NDUFB7 位于模块:", ndufb7_module, "\n")
  
  pdf(paste0("03_results/09_single_cell/hdwgcna/03_hub_network_", ndufb7_module, ".pdf"), width = 8, height = 8)
  p_hub <- HubGeneNetworkPlot(srt, target_modules = ndufb7_module, n_hubs = 10, n_other = 5, edge_alpha = 0.5)
  print(p_hub)
  dev.off()

  module_name <- paste0("CM_OXPHOS-", ndufb7_module)
  if(module_name %in% colnames(srt@meta.data)) {
    pdf(paste0("03_results/09_single_cell/hdwgcna/04_module_vln_", ndufb7_module, ".pdf"), width = 10, height = 6)
    p_vln <- VlnPlot(srt, features = module_name, group.by = "seurat_clusters", split.by = "condition", pt.size = 0)
    print(p_vln)
    dev.off()
  }
} else {
  cat("\n⚠️ NDUFB7 未被分配到模块\n")
}

saveRDS(srt, "03_results/09_single_cell/hdwgcna/srt_hdwgcna_completed.rds")
cat("\n✅ 分析完成！\n")
