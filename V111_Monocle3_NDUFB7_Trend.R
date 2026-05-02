#!/usr/bin/env Rscript
# V111: Monocle3拟时序NDUFB7趋势提取

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

message("========================================")
message("V111: Monocle3 NDUFB7趋势提取")
message("========================================")

# --- 加载cds ---
CDS_FILE <- "03_results/01_seurat_objects/21_monocle3_cds.rds"
if (!file.exists(CDS_FILE)) {
  # 尝试其他路径
  CDS_FILE <- list.files("03_results", pattern = "monocle3.*cds.*\\.rds$", 
                         full.names = TRUE, recursive = TRUE)[1]
}
if (is.na(CDS_FILE) || !file.exists(CDS_FILE)) stop("cds文件不存在")

cds <- readRDS(CDS_FILE)
message("[PASS] CDS加载: ", ncol(cds), " cells × ", nrow(cds), " genes")

# --- 检查NDUFB7 ---
if (!"NDUFB7" %in% rownames(cds)) {
  # 尝试Ensembl ID
  matches <- grep("NDUFB7|99795", rownames(cds), value = TRUE)
  if (length(matches) > 0) {
    target <- matches[1]
    message("[INFO] 使用匹配基因: ", target)
  } else {
    stop("NDUFB7不在cds中")
  }
} else {
  target <- "NDUFB7"
}

# --- 提取NDUFB7沿伪时间表达 ---
message("[EXTRACT] NDUFB7沿伪时间趋势...")

# 获取伪时间
pt <- pseudotime(cds)
pt_clean <- pt[!is.na(pt)]

# 获取NDUFB7表达
expr <- as.numeric(counts(cds)[target, ])
expr_norm <- log1p(expr)

# 拟合趋势（广义加性模型）
df <- data.frame(pseudotime = pt_clean, NDUFB7 = expr_norm[!is.na(pt)])
df <- df[order(df$pseudotime), ]

# 分段统计：早期 vs 晚期
mid_pt <- median(df$pseudotime)
early <- df$NDUFB7[df$pseudotime <= mid_pt]
late <- df$NDUFB7[df$pseudotime > mid_pt]

tt <- t.test(early, late)
cor_test <- cor.test(df$pseudotime, df$NDUFB7, method = "spearman")

message("\n=== Monocle3 统计 ===")
message("伪时间范围: ", round(min(df$pseudotime), 2), " - ", round(max(df$pseudotime), 2))
message("NDUFB7 早期均值: ", round(mean(early), 3))
message("NDUFB7 晚期均值: ", round(mean(late), 3))
message("早期 vs 晚期 t-test p: ", format(tt$p.value, digits = 2, scientific = TRUE))
message("Spearman ρ (伪时间 vs NDUFB7): ", round(cor_test$estimate, 3), 
        " p = ", format(cor_test$p.value, digits = 2, scientific = TRUE))

# 可视化
outdir <- file.path(PROJECT_DIR, "03_results/V111_Monocle3")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

p <- ggplot(df, aes(x = pseudotime, y = NDUFB7)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#31688E") +
  geom_smooth(method = "gam", color = "#FDE725", se = TRUE, linewidth = 1) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("Spearman ρ = ", round(cor_test$estimate, 3),
                         "\np = ", format(cor_test$p.value, digits = 1, scientific = TRUE)),
           hjust = 1.1, vjust = -0.5, size = 3) +
  labs(title = paste0(target, " Expression Along Pseudotime"),
       x = "Pseudotime", y = paste0(target, " Expression (log1p)")) +
  theme_minimal(base_size = 10)

ggsave(file.path(outdir, "V111_NDUFB7_pseudotime.png"), p, width = 6, height = 4, dpi = 300)

# 保存统计
stats <- data.frame(
  N_cells = nrow(df),
  Pseudotime_range = paste0(round(min(df$pseudotime), 2), "-", round(max(df$pseudotime), 2)),
  Early_mean = mean(early),
  Late_mean = mean(late),
  T_test_p = tt$p.value,
  Spearman_rho = cor_test$estimate,
  Spearman_p = cor_test$p.value
)
write.csv(stats, file.path(outdir, "V111_monocle3_stats.csv"), row.names = FALSE)

verdict <- ifelse(cor_test$p.value < 0.05, "PASS", "PARTIAL")
message("\n[", verdict, "] Monocle3趋势提取完成")
message("[DONE] 结果: ", outdir)
