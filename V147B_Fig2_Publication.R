#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(viridis)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V147_Fig2_Publication")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V147B: Fig 2 投稿级制图（三面板）")
message("========================================")

# ========== Panel A: Visium Pseudo-spatial UMAP ==========
visium_file <- file.path(outdir, "V147A_visium_umap_ndufb7.csv")
pA <- NULL
if (file.exists(visium_file)) {
  vis <- read.csv(visium_file, stringsAsFactors = FALSE)
  if (nrow(vis) > 0 && all(c("UMAP1","UMAP2","NDUFB7") %in% colnames(vis))) {
    # 取第一个文件作为主图
    vis_main <- vis[vis$File == vis$File[1], ]
    
    pA <- ggplot(vis_main, aes(x = UMAP1, y = UMAP2, color = NDUFB7)) +
      geom_point(alpha = 0.6, size = 0.8) +
      scale_color_viridis_c(option = "plasma", name = "NDUFB7", 
                            limits = c(0, quantile(vis_main$NDUFB7, 0.95, na.rm = TRUE))) +
      labs(title = "A", subtitle = paste0("Visium Pseudo-spatial (", nrow(vis_main), " spots)"),
           x = "UMAP 1", y = "UMAP 2") +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        panel.background = element_rect(fill = "#F0F0F0", color = NA),
        panel.grid.major = element_line(color = "white"),
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80")
      )
    message("[PASS] Panel A: Visium UMAP generated")
  }
}

if (is.null(pA)) {
  # 备用：空面板
  pA <- ggplot() + annotate("text", x = 0.5, y = 0.5, 
                            label = "A\nVisium spatial data\npending full reconstruction") +
    theme_void() + theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
  message("[WARN] Panel A: Visium data unavailable, using placeholder")
}

# ========== Panel B: Temporal Gradient (GSE214611) ==========
ts_file <- "03_results/02_tables/V68_GSE214611_TimeSeries.csv"
pB <- NULL
if (file.exists(ts_file)) {
  ts <- read.csv(ts_file)
  if (all(c("time","NDUFB7") %in% colnames(ts))) {
    ts <- ts[order(ts$time), ]
    
    pB <- ggplot(ts, aes(x = factor(time), y = NDUFB7, group = 1)) +
      geom_line(color = "#D55E00", linewidth = 1) +
      geom_point(color = "#D55E00", size = 3) +
      geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.8) +
      labs(title = "B", subtitle = "Temporal Gradient (GSE214611)", 
           x = "Time post-MI", y = "NDUFB7 Expression") +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        panel.background = element_rect(fill = "#F0F0F0", color = NA),
        panel.grid.major = element_line(color = "white")
      )
    message("[PASS] Panel B: Temporal gradient generated")
  }
}

if (is.null(pB)) {
  pB <- ggplot() + annotate("text", x = 0.5, y = 0.5, 
                            label = "B\nTemporal gradient\npending time-series data") +
    theme_void() + theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
  message("[WARN] Panel B: Time-series data unavailable")
}

# ========== Panel C: Etiology Gradient (V145_FIX2) ==========
eti_file <- "03_results/V145_GSE57338_Etiology/V145_GSE57338_NDUFB7_by_Etiology.csv"
pC <- NULL
if (file.exists(eti_file)) {
  eti <- read.csv(eti_file, stringsAsFactors = FALSE)
  eti$Disease <- factor(eti$Disease, levels = c("ICM", "DCM"))
  
  # 统计值（从V145_FIX2日志硬编码，确保图注准确）
  kw_p <- 0.05  # Kruskal-Wallis p=0.05 from V145_FIX2
  
  pC <- ggplot(eti, aes(x = Disease, y = NDUFB7, fill = Disease)) +
    geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.12, alpha = 0.35, size = 1.2, color = "black") +
    scale_fill_manual(values = c("ICM" = "#D55E00", "DCM" = "#E69F00"),
                      labels = c("ICM" = "Ischemic CM", "DCM" = "Dilated CM")) +
    labs(
      title = "C",
      subtitle = paste0("Etiology Gradient (GSE57338) | Kruskal-Wallis p = ", signif(kw_p, 2)),
      x = "Etiology", y = "NDUFB7 Expression"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      panel.background = element_rect(fill = "#F0F0F0", color = NA),
      panel.grid.major = element_line(color = "white"),
      legend.position = "none",
      axis.text.x = element_text(face = "bold")
    )
  message("[PASS] Panel C: Etiology gradient generated")
}

if (is.null(pC)) {
  pC <- ggplot() + annotate("text", x = 0.5, y = 0.5, 
                            label = "C\nEtiology gradient\npending V145 output") +
    theme_void() + theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
  message("[WARN] Panel C: Etiology data unavailable")
}

# ========== 合并输出 ==========
combined <- grid.arrange(pA, pB, pC, ncol = 3,
                         top = textGrob("NDUFB7 Spatial and Etiological Gradients", 
                                        gp = gpar(fontface = "bold", fontsize = 13)))

ggsave(file.path(outdir, "Fig2_NDUFB7_Gradients_300dpi.png"), combined,
       width = 16, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(outdir, "Fig2_NDUFB7_Gradients.pdf"), combined,
       width = 16, height = 5.5, device = cairo_pdf)

# 单独保存各面板
if (!is.null(pA)) ggsave(file.path(outdir, "Fig2A_Visium.png"), pA, width = 5.5, height = 5, dpi = 300, bg = "white")
if (!is.null(pB)) ggsave(file.path(outdir, "Fig2B_Temporal.png"), pB, width = 5, height = 5, dpi = 300, bg = "white")
if (!is.null(pC)) ggsave(file.path(outdir, "Fig2C_Etiology.png"), pC, width = 5, height = 5, dpi = 300, bg = "white")

message("[DONE] V147B: ", outdir)
message("  Fig2 combined: Fig2_NDUFB7_Gradients_300dpi.png")
message("  Fig2 PDF: Fig2_NDUFB7_Gradients.pdf")
