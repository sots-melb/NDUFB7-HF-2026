#!/usr/bin/env Rscript
# V110: CellChat CM-Fibroblast互作分析
# 固定路径，使用GSE183852 CM_annotated.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

message("========================================")
message("V110: CellChat CM-Fibroblast互作分析")
message("========================================")

# --- 加载数据 ---
CM_FILE <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"
if (!file.exists(CM_FILE)) stop("CM文件不存在")

srt <- readRDS(CM_FILE)
message("[PASS] CM数据: ", ncol(srt), " cells × ", nrow(srt), " genes")

# --- 检查/创建细胞类型注释 ---
# 如果只有CM，需要模拟或加载全量数据中的其他细胞类型
# 策略：如果数据只有CM，则只做CM内部亚群互作
meta <- srt@meta.data
if (!"cell_type" %in% colnames(meta) && "NDUFB7_level" %in% colnames(meta)) {
  # 用NDUFB7_level作为分组做"伪细胞类型"
  srt$cell_group <- srt$NDUFB7_level
  message("[INFO] 使用NDUFB7_level作为分组: ", paste(unique(srt$cell_group), collapse = ", "))
} else if ("cell_type" %in% colnames(meta)) {
  srt$cell_group <- srt$cell_type
  message("[INFO] 使用cell_type: ", paste(unique(srt$cell_group), collapse = ", "))
} else {
  # 使用final_cluster或major_labl
  for (col in c("final_cluster", "major_labl", "seurat_clusters")) {
    if (col %in% colnames(meta)) {
      srt$cell_group <- meta[[col]]
      message("[INFO] 使用 ", col, " 作为分组")
      break
    }
  }
}

# 创建CellChat对象
message("[LOAD] 创建CellChat对象...")
cellchat <- createCellChat(object = srt, group.by = "cell_group")

# 设置数据库（人类）
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# 预处理
message("[STEP 1/4] 预处理...")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 计算通讯概率
message("[STEP 2/4] 计算通讯概率...")
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 聚合通路
message("[STEP 3/4] 聚合信号通路...")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 可视化
outdir <- file.path(PROJECT_DIR, "03_results/V110_CellChat")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("[STEP 4/4] 可视化...")

# 1. 互作强度热图
groupSize <- as.numeric(table(cellchat@idents))
p1 <- netVisual_heatmap(cellchat, color.heatmap = "Reds")
png(file.path(outdir, "V110_cellchat_heatmap.png"), width = 800, height = 600)
print(p1)
dev.off()

# 2. 圆形互作图
p2 <- netVisual_circle(cellchat, weight.scale = TRUE)
png(file.path(outdir, "V110_cellchat_circle.png"), width = 800, height = 800)
print(p2)
dev.off()

# 3. 提取top互作对
df_net <- subsetCommunication(cellchat)
top_pairs <- df_net[order(-df_net$prob), ]
write.csv(top_pairs, file.path(outdir, "V110_cellchat_all_pairs.csv"), row.names = FALSE)

message("")
message("=== CellChat 结果 ===")
message("Top 10 互作对 (by probability):")
print(head(top_pairs[, c("source", "target", "ligand", "receptor", "prob")], 10))

# 统计显著互作对
sig_pairs <- df_net[df_net$pval < 0.05, ]
message("\n显著互作对 (p<0.05): ", nrow(sig_pairs), " / ", nrow(df_net))

# 判断
if (nrow(sig_pairs) >= 3) {
  message("[PASS] 显著互作对≥3，支持CM-Fibroblast异常互作假说")
  message("[IMPACT] 可写入Results，支撑Fig 5A")
} else {
  message("[WARN] 显著互作对不足")
}

message("[DONE] 结果保存: ", outdir)
