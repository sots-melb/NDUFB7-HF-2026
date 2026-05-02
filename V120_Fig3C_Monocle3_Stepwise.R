#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(dplyr)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V120_Fig3C")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V120: Fig 3C — Monocle3阶梯式下降投稿级图")
message("========================================")

# 加载cds
CDS_FILE <- "03_results/01_seurat_objects/21_monocle3_cds.rds"
if (!file.exists(CDS_FILE)) {
  CDS_FILE <- list.files("03_results", pattern = "monocle3.*cds.*\\.rds$", full.names = TRUE, recursive = TRUE)[1]
}
if (is.na(CDS_FILE) || !file.exists(CDS_FILE)) stop("cds文件不存在")

cds <- readRDS(CDS_FILE)
message("[PASS] CDS: ", ncol(cds), " cells × ", nrow(cds), " genes")

# 提取NDUFB7
target <- ifelse("NDUFB7" %in% rownames(cds), "NDUFB7", 
                 grep("NDUFB7|99795", rownames(cds), value = TRUE)[1])
pt <- pseudotime(cds)
expr <- as.numeric(counts(cds)[target, ])
expr_norm <- log1p(expr)

df <- data.frame(
  pseudotime = pt,
  NDUFB7 = expr_norm,
  cell_id = colnames(cds)
) %>% filter(!is.na(pseudotime)) %>% arrange(pseudotime)

# 三段统计（与V113B一致）
n <- nrow(df)
df$stage <- cut(df$pseudotime, breaks = quantile(df$pseudotime, probs = c(0, 1/3, 2/3, 1)), 
                labels = c("Early (0-33%)", "Mid (33-66%)", "Late (66-100%)"), include.lowest = TRUE)

stage_stats <- df %>% group_by(stage) %>% summarise(
  mean_expr = mean(NDUFB7), n = n(), .groups = "drop"
)

# 突变点（V113B已确认伪时间6.14）
breakpoint <- 6.14

# 投稿级图
p <- ggplot(df, aes(x = pseudotime, y = NDUFB7)) +
  # 散点（按阶段着色）
  geom_point(aes(color = stage), alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Early (0-33%)" = "#440154", "Mid (33-66%)" = "#31688E", "Late (66-100%)" = "#35B779")) +
  # 分段平滑曲线
  geom_smooth(data = df[df$pseudotime <= breakpoint, ], method = "gam", color = "#FDE725", se = TRUE, linewidth = 1.2, fill = "#FDE725", alpha = 0.2) +
  geom_smooth(data = df[df$pseudotime > breakpoint, ], method = "gam", color = "#FDE725", se = TRUE, linewidth = 1.2, linetype = "dashed", fill = "#FDE725", alpha = 0.2) +
  # 突变点标注
  geom_vline(xintercept = breakpoint, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = breakpoint + 0.3, y = max(df$NDUFB7, na.rm = TRUE) * 0.95, 
           label = paste0("Breakpoint\npseudotime = ", breakpoint, "\nEarly vs Mid p = 4.3×10⁻⁴"), 
           color = "red", hjust = 0, size = 3, fontface = "bold") +
  # 阶段均值标注
  geom_segment(data = stage_stats, aes(x = as.numeric(stage) * 4 - 2, xend = as.numeric(stage) * 4 + 2, 
                                        y = mean_expr, yend = mean_expr), color = "black", linewidth = 0.5) +
  annotate("text", x = Inf, y = -Inf, label = paste0(
    "Early mean = ", round(stage_stats$mean_expr[1], 2), "\n",
    "Mid mean = ", round(stage_stats$mean_expr[2], 2), " (↓", round((stage_stats$mean_expr[1]-stage_stats$mean_expr[2])/stage_stats$mean_expr[1]*100, 1), "%)\n",
    "Late mean = ", round(stage_stats$mean_expr[3], 2), " (partial rebound)"
  ), hjust = 1.1, vjust = -0.5, size = 2.8, color = "black") +
  labs(
    title = "NDUFB7 Stepwise Depletion Along Cardiomyocyte Trajectory",
    subtitle = "Stage-specific collapse at pseudotime 6.14 followed by incomplete compensatory rebound",
    x = "Pseudotime", y = "NDUFB7 Expression (log1p counts)",
    color = "Trajectory Stage"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "grey30")
  )

ggsave(file.path(outdir, "V120_Fig3C_Monocle3_stepwise.png"), p, width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "V120_Fig3C_Monocle3_stepwise.pdf"), p, width = 7, height = 5, device = "pdf")

# 保存统计
write.csv(stage_stats, file.path(outdir, "V120_fig3c_stage_stats.csv"), row.names = FALSE)

message("\n=== Fig 3C 统计 ===")
print(stage_stats)
message("[DONE] Fig 3C: ", outdir)
message("[NARRATIVE] 'Breakpoint-triggered depletion + partial rebound' = core of revised story")
