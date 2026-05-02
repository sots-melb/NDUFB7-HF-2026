#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V128_CellChat_Directionality")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V128: CellChat 方向性统计检验")
message("========================================")

# 加载互作数据
chat_file <- "03_results/V113A_CellChat/V113A_cellchat_all_pairs.csv"
if (!file.exists(chat_file)) stop("CellChat结果不存在")
df <- read.csv(chat_file)

# 加载cluster NDUFB7水平
cluster_file <- "03_results/V121_CellChat_Interpretation/V121_cluster_ndufb7_levels.csv"
if (!file.exists(cluster_file)) stop("Cluster NDUFB7水平文件不存在")
clust <- read.csv(cluster_file)

# 映射（V113A使用C编号 = cluster+1）
df$source_cluster <- as.numeric(gsub("C", "", df$source)) - 1
df$target_cluster <- as.numeric(gsub("C", "", df$target)) - 1

# 合并NDUFB7水平
cluster_map <- setNames(clust$NDUFB7_median, clust$seurat_clusters)
df$source_NDUFB7 <- cluster_map[as.character(df$source_cluster)]
df$target_NDUFB7 <- cluster_map[as.character(df$target_cluster)]

# 移除NA（C8等不匹配编号）
df_clean <- df %>% filter(!is.na(source_NDUFB7) & !is.na(target_NDUFB7))

# 方向分类
df_clean$direction <- ifelse(df_clean$source_NDUFB7 < df_clean$target_NDUFB7, "Low_to_High",
                             ifelse(df_clean$source_NDUFB7 > df_clean$target_NDUFB7, "High_to_Low", "Same"))

# --- [1/2] 全局方向性检验 ---
message("\n>>> [1/2] 全局方向性检验")

direction_test <- df_clean %>% 
  group_by(direction) %>% 
  summarise(
    n_pairs = n(),
    mean_prob = mean(prob),
    median_prob = median(prob),
    sd_prob = sd(prob),
    .groups = "drop"
  )
write.csv(direction_test, file.path(outdir, "V128_global_directionality.csv"), row.names = FALSE)
print(direction_test)

# Mann-Whitney: Low→High vs High→Low
low_high <- df_clean$prob[df_clean$direction == "Low_to_High"]
high_low <- df_clean$prob[df_clean$direction == "High_to_Low"]

if (length(low_high) > 5 && length(high_low) > 5) {
  mw_test <- wilcox.test(low_high, high_low, alternative = "greater")
  message("\nMann-Whitney U (Low→High > High→Low): W=", mw_test$statistic, 
          ", p=", format(mw_test$p.value, digits = 2, scientific = TRUE))
  
  if (mw_test$p.value < 0.05) {
    message("[PASS] Low→High interactions are significantly stronger than High→Low")
  } else {
    message("[WARN] No significant directionality asymmetry")
  }
} else {
  message("[SKIP] Insufficient pairs for directionality test")
}

# --- [2/2] ANXA1-FPR1特异性检验 ---
message("\n>>> [2/2] ANXA1-FPR1特异性检验")

# ANXA1-FPR1 vs 其他配体-受体对的概率比较
anxa1_fpr1 <- df_clean %>% filter(ligand == "ANXA1" & receptor == "FPR1")
other_pairs <- df_clean %>% filter(!(ligand == "ANXA1" & receptor == "FPR1"))

if (nrow(anxa1_fpr1) > 0 && nrow(other_pairs) > 0) {
  anxa1_prob <- anxa1_fpr1$prob
  other_prob <- other_pairs$prob
  
  mw_anxa1 <- wilcox.test(anxa1_prob, other_prob, alternative = "greater")
  message("ANXA1-FPR1 median prob: ", round(median(anxa1_prob), 4))
  message("Other pairs median prob: ", round(median(other_prob), 4))
  message("Mann-Whitney (ANXA1-FPR1 > others): p=", format(mw_anxa1$p.value, digits=2, scientific=TRUE))
  
  if (mw_anxa1$p.value < 0.05) {
    message("[PASS] ANXA1-FPR1 is significantly stronger than average LR pair")
  } else {
    message("[PARTIAL] ANXA1-FPR1 not significantly above background")
  }
  
  # 保存
  write.csv(data.frame(
    comparison = "ANXA1-FPR1 vs Others",
    median_anxa1 = median(anxa1_prob),
    median_other = median(other_prob),
    mw_p = mw_anxa1$p.value,
    n_anxa1 = length(anxa1_prob),
    n_other = length(other_prob)
  ), file.path(outdir, "V128_anxa1_specificity.csv"), row.names = FALSE)
}

# 可视化
p <- ggplot(df_clean, aes(x = direction, y = prob, fill = direction)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  scale_fill_manual(values = c("Low_to_High" = "#FDE725", "High_to_Low" = "#440154", "Same" = "#31688E")) +
  labs(title = "Communication Probability by NDUFB7 Directionality",
       subtitle = paste0("Low→High median=", round(median(low_high, na.rm=TRUE), 4), 
                       " | High→Low median=", round(median(high_low, na.rm=TRUE), 4)),
       x = "Signal Direction", y = "Probability") +
  theme_minimal(base_size = 10)

ggsave(file.path(outdir, "V128_directionality_boxplot.png"), p, width = 6, height = 5, dpi = 300)

message("\n[DONE] V128: ", outdir)
