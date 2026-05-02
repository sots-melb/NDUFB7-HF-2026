#!/usr/bin/env Rscript
# V114: 修复CellChat提取 + GSE183852 NDUFB7双峰分布检验（阈值特征）

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(mixtools)  # 双峰检验
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

message("========================================")
message("V114: CellChat修复 + 双峰阈值检验")
message("========================================")

# ========================================
# MODULE 1: CellChat互作矩阵提取（从已计算对象）
# ========================================
message("")
message(">>> [M1] CellChat互作矩阵提取")

CM_FILE <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"
srt <- readRDS(CM_FILE)

# 重新计算（快速，因为数据小）
srt$cell_group <- paste0("C", as.numeric(srt$seurat_clusters) + 1)
cellchat <- createCellChat(object = srt, group.by = "cell_group")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 提取互作矩阵（不画图，避免S4报错）
outdir <- file.path(PROJECT_DIR, "03_results/V114_CellChat")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 提取所有通路
df_net <- subsetCommunication(cellchat)
write.csv(df_net, file.path(outdir, "V114_cellchat_all_pairs.csv"), row.names = FALSE)

# Top统计
top10 <- df_net[order(-df_net$prob), ]
sig <- df_net[df_net$pval < 0.05, ]

message("Top 10 互作对:")
print(head(top10[, c("source","target","ligand","receptor","prob","pval")], 10))
message("\n显著互作对 (p<0.05): ", nrow(sig), " / ", nrow(df_net))

# 提取NDUFB7_high/low相关互作（如果分组可用）
if ("NDUFB7_level" %in% colnames(srt@meta.data)) {
  # 重新用NDUFB7_level做CellChat
  srt$ndufb7_grp <- srt$NDUFB7_level
  cc2 <- createCellChat(object = srt, group.by = "ndufb7_grp")
  cc2@DB <- CellChatDB.human
  cc2 <- subsetData(cc2)
  cc2 <- identifyOverExpressedGenes(cc2)
  cc2 <- identifyOverExpressedInteractions(cc2)
  cc2 <- computeCommunProb(cc2, type = "triMean")
  cc2 <- filterCommunication(cc2, min.cells = 10)
  cc2 <- aggregateNet(cc2)
  
  df2 <- subsetCommunication(cc2)
  write.csv(df2, file.path(outdir, "V114_cellchat_NDUFB7_group_pairs.csv"), row.names = FALSE)
  message("\nNDUFB7分组互作对: ", nrow(df2))
}

# ========================================
# MODULE 2: NDUFB7双峰分布检验（阈值特征）
# ========================================
message("")
message(">>> [M2] GSE183852 NDUFB7双峰分布检验")

ndufb7_expr <- as.numeric(FetchData(srt, vars = "NDUFB7")$NDUFB7)

# 密度图
p1 <- ggplot(data.frame(NDUFB7 = ndufb7_expr), aes(x = NDUFB7)) +
  geom_density(fill = "#31688E", alpha = 0.5) +
  geom_rug(alpha = 0.3) +
  labs(title = "NDUFB7 Expression Distribution in CM (GSE183852)",
       subtitle = paste0("N=", length(ndufb7_expr), " | Testing for bimodal pattern"),
       x = "NDUFB7 Expression", y = "Density") +
  theme_minimal(base_size = 10)

ggsave(file.path(outdir, "V114_NDUFB7_density.png"), p1, width = 6, height = 4, dpi = 300)

# 双峰检验（如果mixtools可用）
bimodal_test <- tryCatch({
  mix <- normalmixEM(ndufb7_expr, k = 2, maxit = 100, epsilon = 1e-6)
  # 计算双峰系数 (BC): (μ1-μ2)^2 / (σ1^2+σ2^2)
  bc <- (mix$mu[1] - mix$mu[2])^2 / (mix$sigma[1]^2 + mix$sigma[2]^2)
  list(mu = mix$mu, sigma = mix$sigma, lambda = mix$lambda, bc = bc, loglik = mix$loglik)
}, error = function(e) {
  message("[WARN] mixtools失败: ", conditionMessage(e))
  NULL
})

if (!is.null(bimodal_test)) {
  message("\n=== 双峰混合模型 ===")
  message("组分1: μ=", round(bimodal_test$mu[1], 3), 
          " σ=", round(bimodal_test$sigma[1], 3),
          " λ=", round(bimodal_test$lambda[1], 3))
  message("组分2: μ=", round(bimodal_test$mu[2], 3), 
          " σ=", round(bimodal_test$sigma[2], 3),
          " λ=", round(bimodal_test$lambda[2], 3))
  message("双峰系数 (BC): ", round(bimodal_test$bc, 3))
  message("[INTERPRETATION] BC>2 强烈支持双峰/阈值特征")
  
  # 叠加拟合曲线
  x_seq <- seq(min(ndufb7_expr), max(ndufb7_expr), length.out = 200)
  y1 <- dnorm(x_seq, bimodal_test$mu[1], bimodal_test$sigma[1]) * bimodal_test$lambda[1]
  y2 <- dnorm(x_seq, bimodal_test$mu[2], bimodal_test$sigma[2]) * bimodal_test$lambda[2]
  
  fit_df <- data.frame(x = x_seq, y1 = y1, y2 = y2, total = y1 + y2)
  
  p2 <- p1 + 
    geom_line(data = fit_df, aes(x = x, y = total), color = "red", linewidth = 1) +
    geom_line(data = fit_df, aes(x = x, y = y1), color = "blue", linetype = "dashed") +
    geom_line(data = fit_df, aes(x = x, y = y2), color = "green", linetype = "dashed") +
    annotate("text", x = Inf, y = Inf,
             label = paste0("BC=", round(bimodal_test$bc, 2)),
             hjust = 1.1, vjust = 1.5, size = 3, color = "red")
  
  ggsave(file.path(outdir, "V114_NDUFB7_bimodal_fit.png"), p2, width = 6, height = 4, dpi = 300)
} else {
  # 简单直方图检验
  h <- hist(ndufb7_expr, breaks = 30, plot = FALSE)
  # 找局部最小值（valley）
  valleys <- which(diff(sign(diff(h$density))) == 2) + 1
  message("\n[MANUAL] 密度 valleys 位置: ", paste(round(h$mids[valleys], 3), collapse = ", "))
  message("[INTERPRETATION] 如果存在明显valley，支持双峰/阈值")
}

# 分位数检验（寻找断崖点）
quants <- quantile(ndufb7_expr, probs = seq(0.1, 0.9, by = 0.1))
message("\n=== 分位数检验 ===")
print(round(quants, 3))

# 找最大密度下降区间
dens <- density(ndufb7_expr)
max_drop_idx <- which.min(diff(dens$y))
threshold_candidate <- dens$x[max_drop_idx]
message("最大密度下降点 (阈值候选): ", round(threshold_candidate, 3))

# 保存
stats <- data.frame(
  Metric = c("N_cells", "Mean", "Median", "SD", "CV", "BC", "Threshold_candidate"),
  Value = c(length(ndufb7_expr), mean(ndufb7_expr), median(ndufb7_expr),
            sd(ndufb7_expr), sd(ndufb7_expr)/mean(ndufb7_expr),
            ifelse(is.null(bimodal_test), NA, bimodal_test$bc),
            threshold_candidate)
)
write.csv(stats, file.path(outdir, "V114_NDUFB7_distribution_stats.csv"), row.names = FALSE)

message("\n[DONE] V114完成: ", outdir)
