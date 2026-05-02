#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V144_Fig3_Publication")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V144: Fig 3 投稿级制图")
message("========================================")

# 读取V141数据
df <- read.csv("03_results/V141_GSE57338_Bulk_Fix/V141_bulk_ferroptosis.csv", stringsAsFactors = FALSE)
df$NDUFB7_group <- factor(df$NDUFB7_group, levels = c("NDUFB7_High", "NDUFB7_Low"))

# 统计值（从V141日志提取，硬编码确保图注准确）
rho_val <- -0.229
p_spear <- 4.6e-05
p_wilcox <- 3.9e-03
rho_acsl4 <- -0.157
p_acsl4 <- 5.6e-03

# 配色方案（Nature风格）
col_high <- "#0072B2"  # 蓝
col_low  <- "#D55E00"  # 橙
col_bg   <- "#F0F0F0"  # 浅灰背景

# ========== Panel A: 散点图 + Spearman ==========
pA <- ggplot(df, aes(x = NDUFB7, y = ferroptosis_score)) +
  geom_point(aes(color = NDUFB7_group), alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("NDUFB7_High" = col_high, "NDUFB7_Low" = col_low),
                     labels = c("NDUFB7_High" = "High", "NDUFB7_Low" = "Low")) +
  labs(
    title = "A",
    x = "NDUFB7 Expression",
    y = "Ferroptosis Score (z-mean)",
    color = "NDUFB7"
  ) +
  annotate("text", x = Inf, y = -Inf, 
           label = paste0("Spearman ρ = ", rho_val, "\np = ", format(p_spear, digits=1, scientific=TRUE)),
           hjust = 1.1, vjust = -0.5, size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    panel.background = element_rect(fill = col_bg, color = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "gray80")
  )

# ========== Panel B: Boxplot + Wilcoxon ==========
pB <- ggplot(df, aes(x = NDUFB7_group, y = ferroptosis_score, fill = NDUFB7_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5, color = "black") +
  scale_fill_manual(values = c("NDUFB7_High" = col_high, "NDUFB7_Low" = col_low),
                    labels = c("NDUFB7_High" = "High", "NDUFB7_Low" = "Low")) +
  labs(
    title = "B",
    x = "NDUFB7 Group (median split)",
    y = "Ferroptosis Score",
    fill = "NDUFB7"
  ) +
  annotate("text", x = 1.5, y = Inf,
           label = paste0("Wilcoxon p = ", format(p_wilcox, digits=1, scientific=TRUE)),
           vjust = 1.5, size = 3.5, fontface = "italic") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    panel.background = element_rect(fill = col_bg, color = NA),
    panel.grid.major = element_line(color = "white"),
    legend.position = "none"
  )

# ========== Panel C: ACSL4/GPX4 Ratio ==========
has_ratio <- "ACSL4_GPX4_ratio" %in% names(df) && sum(is.finite(df$ACSL4_GPX4_ratio)) > 10

if (has_ratio) {
  pC <- ggplot(df, aes(x = NDUFB7, y = ACSL4_GPX4_ratio)) +
    geom_point(aes(color = NDUFB7_group), alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dotdash", linewidth = 0.8) +
    scale_color_manual(values = c("NDUFB7_High" = col_high, "NDUFB7_Low" = col_low)) +
    labs(
      title = "C",
      x = "NDUFB7 Expression",
      y = "ACSL4 / GPX4 Ratio",
      color = "NDUFB7"
    ) +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("Spearman ρ = ", rho_acsl4, "\np = ", format(p_acsl4, digits=1, scientific=TRUE)),
             hjust = 1.1, vjust = -0.5, size = 3.5, fontface = "italic") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      panel.background = element_rect(fill = col_bg, color = NA),
      panel.grid.major = element_line(color = "white"),
      legend.position = c(0.85, 0.85),
      legend.background = element_rect(fill = "white", color = "gray80")
    )
} else {
  # 如果没有ratio数据，画一个空面板占位
  pC <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, label = "C\nACSL4/GPX4 data pending") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
}

# ========== 合并输出 ==========
combined <- grid.arrange(pA, pB, pC, ncol = 3, 
                         top = textGrob("NDUFB7 Depletion Activates Ferroptosis (GSE57338, n=313)", 
                                        gp = gpar(fontface = "bold", fontsize = 14)))

ggsave(file.path(outdir, "Fig3_NDUFB7_Ferroptosis_300dpi.png"), combined, 
       width = 16, height = 6, dpi = 300, bg = "white")
ggsave(file.path(outdir, "Fig3_NDUFB7_Ferroptosis.pdf"), combined, 
       width = 16, height = 6, device = cairo_pdf)

# 单独保存各面板（方便排版调整）
ggsave(file.path(outdir, "Fig3A_scatter.png"), pA, width = 5.5, height = 5, dpi = 300, bg = "white")
ggsave(file.path(outdir, "Fig3B_boxplot.png"), pB, width = 5, height = 5, dpi = 300, bg = "white")
if (has_ratio) ggsave(file.path(outdir, "Fig3C_ratio.png"), pC, width = 5.5, height = 5, dpi = 300, bg = "white")

message("[DONE] V144: ", outdir)
message("  Fig3 combined: Fig3_NDUFB7_Ferroptosis_300dpi.png")
message("  Fig3 PDF: Fig3_NDUFB7_Ferroptosis.pdf")
