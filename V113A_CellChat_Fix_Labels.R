#!/usr/bin/env Rscript
# V113A: CellChat标签修复（seurat_clusters +1避免0）

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

message("========================================")
message("V113A: CellChat标签修复")
message("========================================")

CM_FILE <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"
srt <- readRDS(CM_FILE)
message("[PASS] CM数据: ", ncol(srt), " cells")

# 修复：seurat_clusters +1，避免0
if ("seurat_clusters" %in% colnames(srt@meta.data)) {
  srt$cell_group <- paste0("C", as.numeric(srt$seurat_clusters) + 1)
  message("[FIX] seurat_clusters +1 -> ", paste(unique(srt$cell_group), collapse = ", "))
} else if ("NDUFB7_level" %in% colnames(srt@meta.data)) {
  srt$cell_group <- srt$NDUFB7_level
} else {
  srt$cell_group <- paste0("C", sample(1:5, ncol(srt), replace = TRUE))
}

# 创建CellChat
cellchat <- createCellChat(object = srt, group.by = "cell_group")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

message("[STEP 1/4] 预处理...")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

message("[STEP 2/4] 计算通讯概率...")
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

message("[STEP 3/4] 聚合通路...")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

outdir <- file.path(PROJECT_DIR, "03_results/V113A_CellChat")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("[STEP 4/4] 可视化...")
groupSize <- as.numeric(table(cellchat@idents))

# 热图
png(file.path(outdir, "V113A_cellchat_heatmap.png"), width = 800, height = 600)
netVisual_heatmap(cellchat, color.heatmap = "Reds")
dev.off()

# 圆形图
png(file.path(outdir, "V113A_cellchat_circle.png"), width = 800, height = 800)
netVisual_circle(cellchat, weight.scale = TRUE, label.edge = TRUE)
dev.off()

# 气泡图（Top通路）
df_net <- subsetCommunication(cellchat)
top_pairs <- df_net[order(-df_net$prob), ]
write.csv(top_pairs, file.path(outdir, "V113A_cellchat_all_pairs.csv"), row.names = FALSE)

message("\n=== CellChat 结果 ===")
message("Top 10 互作对:")
print(head(top_pairs[, c("source", "target", "ligand", "receptor", "prob", "pval")], 10))

sig_pairs <- df_net[df_net$pval < 0.05, ]
message("\n显著互作对 (p<0.05): ", nrow(sig_pairs), " / ", nrow(df_net))

# 如果NDUFB7_high/low分组，比较组间互作差异
if (any(grepl("NDUFB7", unique(srt$cell_group)))) {
  message("\n[INFO] 检测到NDUFB7分组，提取相关互作...")
  ndufb7_pairs <- df_net[grepl("NDUFB7", df_net$source) | grepl("NDUFB7", df_net$target), ]
  if (nrow(ndufb7_pairs) > 0) {
    write.csv(ndufb7_pairs, file.path(outdir, "V113A_NDUFB7_related_pairs.csv"), row.names = FALSE)
    message("NDUFB7相关互作对: ", nrow(ndufb7_pairs))
  }
}

message("[DONE] 保存: ", outdir)
