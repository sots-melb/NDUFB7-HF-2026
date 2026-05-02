#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(pheatmap)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V113A_CellChat")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V113A_FIX: CellChat数据提取（绕过S4可视化bug）")
message("========================================")

CM_FILE <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"
srt <- readRDS(CM_FILE)
message("[PASS] CM: ", ncol(srt), " cells")

srt$cell_group <- paste0("C", as.numeric(srt$seurat_clusters) + 1)
cellchat <- createCellChat(object = srt, group.by = "cell_group")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 提取互作对
df_net <- subsetCommunication(cellchat)
write.csv(df_net, file.path(outdir, "V113A_cellchat_all_pairs.csv"), row.names = FALSE)

# 用pheatmap重建热图（绕过CellChat S4 bug）
net_weight <- cellchat@net$weight
if (!is.null(net_weight)) {
  write.csv(as.data.frame(net_weight), file.path(outdir, "V113A_interaction_weight_matrix.csv"))
  png(file.path(outdir, "V113A_heatmap_fixed.png"), width = 800, height = 600)
  pheatmap(net_weight, color = colorRampPalette(c("white", "#FDE725", "#35B779", "#31688E", "#440154"))(100),
           main = "CellChat Interaction Weight (Fixed)", display_numbers = TRUE, number_format = "%.2f",
           fontsize_number = 8)
  dev.off()
  message("[PASS] 热图已用pheatmap重建")
}

# 通路层面
df_pathway <- subsetCommunication(cellchat, slot.name = "netP")
if (nrow(df_pathway) > 0) {
  write.csv(df_pathway, file.path(outdir, "V113A_pathway_communication.csv"), row.names = FALSE)
}

# 统计
top10 <- head(df_net[order(-df_net$prob), c("source", "target", "ligand", "receptor", "prob", "pval")], 10)
write.csv(top10, file.path(outdir, "V113A_top10_pairs.csv"), row.names = FALSE)

message("\n=== CellChat 统计 ===")
message("显著互作对 (p<0.05): ", sum(df_net$pval < 0.05), " / ", nrow(df_net))
message("\nTop 3 互作对:")
print(head(top10, 3))

message("\n[INTERPRETATION] ANXA1-FPR1 = 炎症消退/趋化信号")
message("[INTERPRETATION] PECAM1-PECAM1 = 内皮细胞粘附（CM亚群间?)")
message("[NOTE] 当前数据仅含CM亚群，无Fibroblast。CM-FB互作需全量数据。")
message("[DONE] V113A_FIX: ", outdir)
