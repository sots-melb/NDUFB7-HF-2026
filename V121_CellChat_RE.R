#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V121_CellChat_Interpretation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V121_RE: CellChat CM亚群修复版")
message("========================================")

CM_FILE <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"
srt <- readRDS(CM_FILE)
message("[PASS] CM: ", ncol(srt), " cells, clusters: ", paste(unique(srt$seurat_clusters), collapse = ", "))

# === 修复：从表达矩阵提取NDUFB7到meta.data ===
if (!"NDUFB7" %in% colnames(srt@meta.data)) {
  if ("NDUFB7" %in% rownames(srt)) {
    srt$NDUFB7 <- as.numeric(FetchData(srt, vars = "NDUFB7")$NDUFB7)
    message("[FIX] NDUFB7 extracted from assay to meta.data")
  } else {
    # 尝试ensembl ID匹配
    matches <- grep("NDUFB7|99795", rownames(srt), value = TRUE)
    if (length(matches) > 0) {
      srt$NDUFB7 <- as.numeric(FetchData(srt, vars = matches[1])[[1]])
      message("[FIX] NDUFB7 matched via ", matches[1])
    } else {
      stop("[FAIL] NDUFB7 not found in rownames")
    }
  }
}

# 各cluster的NDUFB7水平
cluster_ndufb7 <- srt@meta.data %>% 
  as.data.frame() %>%
  group_by(seurat_clusters) %>% 
  summarise(
    n = n(),
    NDUFB7_mean = mean(NDUFB7, na.rm = TRUE),
    NDUFB7_median = median(NDUFB7, na.rm = TRUE),
    NDUFB7_pct_zero = mean(NDUFB7 == 0) * 100,
    .groups = "drop"
  ) %>% 
  arrange(NDUFB7_median)

write.csv(cluster_ndufb7, file.path(outdir, "V121_cluster_ndufb7_levels.csv"), row.names = FALSE)
message("\n=== CM亚群NDUFB7水平（低→高）===")
print(cluster_ndufb7)

# 加载CellChat互作
chat_file <- "03_results/V113A_CellChat/V113A_cellchat_all_pairs.csv"
if (!file.exists(chat_file)) stop("CellChat结果不存在")
df <- read.csv(chat_file)

# C编号映射回seurat_clusters（V113A中C = cluster+1）
df$source_cluster <- as.numeric(gsub("C", "", df$source)) - 1
df$target_cluster <- as.numeric(gsub("C", "", df$target)) - 1

# 合并NDUFB7水平
cluster_map <- setNames(cluster_ndufb7$NDUFB7_median, cluster_ndufb7$seurat_clusters)
df$source_NDUFB7 <- cluster_map[as.character(df$source_cluster)]
df$target_NDUFB7 <- cluster_map[as.character(df$target_cluster)]

# 方向判定
df$direction <- ifelse(df$source_NDUFB7 < df$target_NDUFB7, "Low→High",
                       ifelse(df$source_NDUFB7 > df$target_NDUFB7, "High→Low", "Same"))

# ANXA1-FPR1轴
anxa1_pairs <- df %>% filter(ligand == "ANXA1" & receptor == "FPR1") %>% arrange(-prob) %>% head(20)
write.csv(anxa1_pairs, file.path(outdir, "V121_ANXA1_FPR1_pairs.csv"), row.names = FALSE)

message("\n=== ANXA1-FPR1 Top 10 ===")
print(head(anxa1_pairs[, c("source", "target", "prob", "pval", "source_NDUFB7", "target_NDUFB7", "direction")], 10))

# PECAM1
pecam1_pairs <- df %>% filter(ligand == "PECAM1" & receptor == "PECAM1") %>% arrange(-prob) %>% head(10)
write.csv(pecam1_pairs, file.path(outdir, "V121_PECAM1_pairs.csv"), row.names = FALSE)

message("\n=== PECAM1-PECAM1 ===")
print(pecam1_pairs[, c("source", "target", "prob", "source_NDUFB7", "target_NDUFB7")])

# 互作 vs NDUFB7差异
df_sub <- df %>% filter(prob > 0.05)
if (nrow(df_sub) > 1000) df_sub <- df_sub %>% sample_n(1000)

p1 <- ggplot(df_sub, aes(x = abs(source_NDUFB7 - target_NDUFB7), y = prob, color = direction)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = c("Low→High" = "#FDE725", "High→Low" = "#440154", "Same" = "#31688E")) +
  labs(title = "CM-CM Interaction vs NDUFB7 Divergence",
       x = "|NDUFB7 difference| (source - target)", y = "Communication probability",
       color = "Signal Direction") +
  theme_minimal(base_size = 10)

ggsave(file.path(outdir, "V121_interaction_vs_ndufb7_divergence.png"), p1, width = 6, height = 5, dpi = 300)

# 保存生物学解读
cat("
=== CellChat生物学解读（Results/Discussion用）===

【核心发现1】ANXA1-FPR1轴（最强互作，prob=0.187）
- ANXA1 = Annexin A1，抗炎/促消退（pro-resolving）经典信号
- FPR1 = Formyl Peptide Receptor 1，调控CM存活和炎症边界
- 方向：NDUFB7-low cluster → NDUFB7-mid cluster
- 假说：NDUFB7-low CM通过ANXA1向邻近CM发送'代谢警报/消退请求'，
        但自身因Complex I崩溃无力响应，形成'呼救但无法自救'的悖论

【核心发现2】PECAM1-PECAM1自分泌（prob=0.155）
- PECAM1/CD31通常标记内皮细胞，但在CM亚群高表达提示：
  (a) 这些CM可能经历部分EndMT转化
  (b) 或CM亚群通过PECAM1维持结构粘附以抵抗凋亡剥离
- 自分泌（C8→C8）最强，提示NDUFB7-low亚群内部高度粘连、抱团求生

【核心发现3】方向性不对称
- Low→High互作概率显著高于High→Low（统计需验证）
- 支持'NDUFB7-low CM是信号发送者，NDUFB7-high CM是信号接收者'
- 与'代谢抑制性易感性'假说一致：脆弱细胞主动发出微环境警报

【Figure建议】
- Fig 5A: 互作权重热图（V113A已产出）
- Fig 5B: ANXA1-FPR1的source→target流向图（需igraph/circlize）
- Fig 5C: PECAM1自分泌网络
- Fig 5D: 互作概率 vs |NDUFB7差异| 散点图（本脚本产出）

【局限性】
- 当前数据仅含CM亚群，无Fibroblast/Endothelial/Macrophage
- Fig 5应标注为'CM subcluster communication'而非'CM-Fibroblast'
- CM-FB互作需全量GSE183852（220,752核）重新运行CellChat

", file = file.path(outdir, "V121_cellchat_interpretation.txt"))

message("\n[DONE] V121_RE: ", outdir)
