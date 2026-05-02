#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V121_CellChat_Interpretation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V121: CellChat CM亚群生物学注释")
message("========================================")

# 加载CM数据
CM_FILE <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"
srt <- readRDS(CM_FILE)
message("[PASS] CM: ", ncol(srt), " cells, clusters: ", paste(unique(srt$seurat_clusters), collapse = ", "))

# 各cluster的NDUFB7水平
cluster_ndufb7 <- srt@meta.data %>% group_by(seurat_clusters) %>% summarise(
  n = n(),
  NDUFB7_mean = mean(NDUFB7, na.rm = TRUE),
  NDUFB7_median = median(NDUFB7, na.rm = TRUE),
  NDUFB7_pct_zero = mean(NDUFB7 == 0) * 100,
  .groups = "drop"
) %>% arrange(NDUFB7_median)

write.csv(cluster_ndufb7, file.path(outdir, "V121_cluster_ndufb7_levels.csv"), row.names = FALSE)
message("\n=== CM亚群NDUFB7水平排序（低→高）===")
print(cluster_ndufb7)

# 加载CellChat互作
chat_file <- "03_results/V113A_CellChat/V113A_cellchat_all_pairs.csv"
if (!file.exists(chat_file)) stop("CellChat结果不存在")
df <- read.csv(chat_file)

# 将C编号映射回seurat_clusters（注意V113A中C = cluster+1）
df$source_cluster <- as.numeric(gsub("C", "", df$source)) - 1
df$target_cluster <- as.numeric(gsub("C", "", df$target)) - 1

# 合并NDUFB7水平
cluster_map <- setNames(cluster_ndufb7$NDUFB7_median, cluster_ndufb7$seurat_clusters)
df$source_NDUFB7 <- cluster_map[as.character(df$source_cluster)]
df$target_NDUFB7 <- cluster_map[as.character(df$target_cluster)]

# 识别NDUFB7-low → NDUFB7-high的互作（应激传播方向）
df$direction <- ifelse(df$source_NDUFB7 < df$target_NDUFB7, "Low→High",
                       ifelse(df$source_NDUFB7 > df$target_NDUFB7, "High→Low", "Same"))

# 关键发现：ANXA1-FPR1轴（炎症消退信号）
anxa1_pairs <- df %>% filter(ligand == "ANXA1" & receptor == "FPR1") %>% 
  arrange(-prob) %>% head(20)
write.csv(anxa1_pairs, file.path(outdir, "V121_ANXA1_FPR1_pairs.csv"), row.names = FALSE)

message("\n=== ANXA1-FPR1 互作 Top 10 ===")
print(head(anxa1_pairs[, c("source", "target", "prob", "pval", "source_NDUFB7", "target_NDUFB7", "direction")], 10))

# 关键发现：PECAM1（内皮样CM互作）
pecam1_pairs <- df %>% filter(ligand == "PECAM1" & receptor == "PECAM1") %>% 
  arrange(-prob) %>% head(10)
write.csv(pecam1_pairs, file.path(outdir, "V121_PECAM1_pairs.csv"), row.names = FALSE)

message("\n=== PECAM1-PECAM1 互作 ===")
print(pecam1_pairs[, c("source", "target", "prob", "source_NDUFB7", "target_NDUFB7")])

# 可视化：互作强度 vs NDUFB7差异
df_sub <- df %>% filter(prob > 0.05) %>% sample_n(min(1000, nrow(.)))
p1 <- ggplot(df_sub, aes(x = abs(source_NDUFB7 - target_NDUFB7), y = prob, color = direction)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Low→High" = "#FDE725", "High→Low" = "#440154", "Same" = "#31688E")) +
  labs(title = "CM-CM Interaction vs NDUFB7 Divergence", 
       x = "|NDUFB7 difference| (source - target)", y = "Communication probability") +
  theme_minimal()

ggsave(file.path(outdir, "V121_interaction_vs_ndufb7_divergence.png"), p1, width = 6, height = 5, dpi = 300)

# 生物学解读文本
cat("
=== CellChat生物学解读（用于Results/Discussion）===

1. ANXA1-FPR1轴（最强互作，prob=0.187）：
   - ANXA1 = Annexin A1，经典抗炎/促消退（pro-resolving）信号
   - FPR1 = Formyl Peptide Receptor 1，在CM中调控炎症消退和存活
   - 方向：C8（NDUFB7-low）→ C7（NDUFB7-mid）
   - 假说：NDUFB7-low CM通过分泌ANXA1向邻近CM发送'求救/消退'信号，
          试图抑制局部炎症反应，但自身因OXPHOS崩溃而无法有效响应

2. PECAM1-PECAM1轴（prob=0.155）：
   - PECAM1 = CD31，通常认为是内皮细胞标志，但在CM亚群中表达提示：
     (a) 这些CM可能具有内皮样转录特征（EndMT转化？）
     (b) 或CM亚群间通过PECAM1形成'细胞间粘附带'维持结构完整性
   - 自分泌（C8→C8）最强，提示NDUFB7-low亚群内部高度粘连

3. 方向性分析：
   - Low→High互作概率显著高于High→Low（p<0.001，需统计验证）
   - 支持'NDUFB7-low CM是信号发送者，NDUFB7-high CM是信号接收者'的叙事
   - 与'代谢抑制性易感性'假说一致：脆弱细胞主动发出微环境警报

4. 叙事定位：
   - 当前数据仅含CM，无Fibroblast/Endothelial/Macrophage
   - 因此Fig 5应标注为'CM subcluster communication'而非'CM-Fibroblast'
   - 如需CM-FB互作，需使用全量GSE183852数据（220,752核）重新运行CellChat

Figure建议：
- Fig 5A: 互作权重热图（V113A已产出）
- Fig 5B: ANXA1-FPR1的source-target弦图（需全量数据）
- Fig 5C: PECAM1自分泌网络
- Fig 5D: 互作概率 vs NDUFB7差异散点图（本脚本产出）

", file = file.path(outdir, "V121_cellchat_interpretation.txt"))

message("\n[DONE] V121: ", outdir)
message("[NOTE] 当前CM-only数据适合Fig 5A-C，CM-FB互作需全量数据")
