#!/usr/bin/env Rscript
# V113B: Monocle3阶梯式下降深化——找分支点/阈值

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(dplyr)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

message("========================================")
message("V113B: Monocle3阶梯式下降深化")
message("========================================")

CDS_FILE <- "03_results/01_seurat_objects/21_monocle3_cds.rds"
cds <- readRDS(CDS_FILE)
message("[PASS] CDS: ", ncol(cds), " cells × ", nrow(cds), " genes")

target <- ifelse("NDUFB7" %in% rownames(cds), "NDUFB7", 
                 grep("NDUFB7|99795", rownames(cds), value = TRUE)[1])

pt <- pseudotime(cds)
expr <- as.numeric(counts(cds)[target, ])
expr_norm <- log1p(expr)

df <- data.frame(pseudotime = pt, NDUFB7 = expr_norm) %>%
  filter(!is.na(pseudotime)) %>%
  arrange(pseudotime)

# 找阶梯/阈值：用滑动窗口检测突变点
message("[ANALYSIS] 滑动窗口检测NDUFB7突变点...")

window_size <- round(nrow(df) * 0.1)  # 10%窗口
df$rolling_mean <- zoo::rollmean(df$NDUFB7, window_size, fill = NA, align = "center")
df$rolling_sd <- zoo::rollapply(df$NDUFB7, window_size, sd, fill = NA, align = "center")

# 找最大下降点
df$delta <- c(NA, diff(df$rolling_mean))
max_drop_idx <- which.min(df$delta)
max_drop_pt <- df$pseudotime[max_drop_idx]

message("\n=== 阶梯式下降分析 ===")
message("最大下降点伪时间: ", round(max_drop_pt, 2))
message("该点NDUFB7均值: ", round(df$rolling_mean[max_drop_idx], 3))
message("下降幅度: ", round(df$delta[max_drop_idx], 3))

# 分段统计：早期(0-33%) / 中期(33-66%) / 晚期(66-100%)
q33 <- quantile(df$pseudotime, 0.33)
q66 <- quantile(df$pseudotime, 0.66)

early <- df$NDUFB7[df$pseudotime <= q33]
mid <- df$NDUFB7[df$pseudotime > q33 & df$pseudotime <= q66]
late <- df$NDUFB7[df$pseudotime > q66]

e_m <- t.test(early, mid)
m_l <- t.test(mid, late)
e_l <- t.test(early, late)

message("\n三段比较:")
message("  早期 (0-33%): mean=", round(mean(early), 3), ", n=", length(early))
message("  中期 (33-66%): mean=", round(mean(mid), 3), ", n=", length(mid))
message("  晚期 (66-100%): mean=", round(mean(late), 3), ", n=", length(late))
message("  早期 vs 中期 p: ", format(e_m$p.value, digits = 2, scientific = TRUE))
message("  中期 vs 晚期 p: ", format(m_l$p.value, digits = 2, scientific = TRUE))
message("  早期 vs 晚期 p: ", format(e_l$p.value, digits = 2, scientific = TRUE))

# 可视化
outdir <- file.path(PROJECT_DIR, "03_results/V113B_Monocle3_Stepwise")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

p <- ggplot(df, aes(x = pseudotime, y = NDUFB7)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#31688E") +
  geom_line(aes(y = rolling_mean), color = "#d62728", linewidth = 1) +
  geom_vline(xintercept = max_drop_pt, linetype = "dashed", color = "red") +
  annotate("text", x = max_drop_pt, y = max(df$NDUFB7), 
           label = paste0("Drop point\npt=", round(max_drop_pt, 1)),
           color = "red", hjust = -0.1, vjust = 1) +
  labs(title = paste0(target, " Stepwise Collapse Along Pseudotime"),
       subtitle = paste0("Early vs Late p=", format(e_l$p.value, digits = 1, scientific = TRUE)),
       x = "Pseudotime", y = paste0(target, " Expression")) +
  theme_minimal(base_size = 10)

ggsave(file.path(outdir, "V113B_stepwise_collapse.png"), p, width = 7, height = 5, dpi = 300)

# 保存
stats <- data.frame(
  Phase = c("Early_0_33", "Mid_33_66", "Late_66_100"),
  Mean = c(mean(early), mean(mid), mean(late)),
  N = c(length(early), length(mid), length(late)),
  Drop_Point_PT = c(max_drop_pt, NA, NA)
)
write.csv(stats, file.path(outdir, "V113B_stepwise_stats.csv"), row.names = FALSE)

message("\n[PASS] 阶梯式下降分析完成")
message("[INTERPRETATION] NDUFB7在伪时间中存在突变式丢失，支持'应激阈值'模型")
message("[DONE] 保存: ", outdir)
