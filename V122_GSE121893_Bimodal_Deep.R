#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mclust)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V122_GSE121893")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V122: GSE121893双峰深化（Supplementary）")
message("========================================")

f121893 <- "~/Downloads/GSE121893_human_heart_sc_umi.csv.gz"
if (!file.exists(f121893)) { message("[FAIL] GSE121893 not found"); quit(status=1) }

dt <- fread(f121893, header = TRUE)
gene_col <- colnames(dt)[1]
idx <- which(toupper(dt[[gene_col]]) == "NDUFB7")[1]
if (is.na(idx)) { message("[FAIL] NDUFB7 not found"); quit(status=1) }

expr <- as.numeric(dt[idx, -1, with = FALSE])
expr <- expr[!is.na(expr)]
message("[PASS] GSE121893 NDUFB7: n=", length(expr), " zero_pct=", round(mean(expr==0)*100,1), "%")

# 双峰模型详细参数
mod <- densityMclust(expr, G = 2, verbose = FALSE)
peaks <- sort(mod$parameters$mean)
proportions <- mod$parameters$pro[order(mod$parameters$mean)]

df_detail <- data.frame(
  peak = c("Silent/Low", "High/Retained"),
  mean = peaks,
  sd = sqrt(mod$parameters$variance$sigmasq[order(mod$parameters$mean)]),
  proportion = proportions,
  stringsAsFactors = FALSE
)
write.csv(df_detail, file.path(outdir, "V121_bimodal_parameters.csv"), row.names = FALSE)

message("\n=== GSE121893 双峰参数 ===")
print(df_detail)

# 密度图（投稿级）
df_plot <- data.frame(NDUFB7 = expr)
p <- ggplot(df_plot, aes(x = NDUFB7)) +
  geom_density(fill = "#31688E", alpha = 0.3, color = "#440154", linewidth = 1) +
  geom_vline(xintercept = peaks, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = peaks[1], y = Inf, label = paste0("Peak 1\nμ=", round(peaks[1],2), "\n", round(proportions[1]*100,1), "%"), 
           vjust = 1.5, color = "red", size = 3) +
  annotate("text", x = peaks[2], y = Inf, label = paste0("Peak 2\nμ=", round(peaks[2],2), "\n", round(proportions[2]*100,1), "%"), 
           vjust = 1.5, color = "red", size = 3) +
  labs(
    title = "NDUFB7 Bimodal Distribution in Human Heart scRNA-seq (GSE121893)",
    subtitle = paste0("Mixture model G=2 | All-or-none index = ", round((mean(expr==0) + mean(expr >= quantile(expr[expr>0], 0.9)))*100, 1), "%"),
    x = "NDUFB7 Expression (UMI)", y = "Density"
  ) +
  theme_minimal(base_size = 10)

ggsave(file.path(outdir, "V122_GSE121893_bimodal_density.png"), p, width = 6, height = 4, dpi = 300)
ggsave(file.path(outdir, "V122_GSE121893_bimodal_density.pdf"), p, width = 6, height = 4, device = "pdf")

message("[DONE] V122: ", outdir)
message("[USAGE] 作为Supplementary Figure S3: 'Bimodal distribution in independent cohort'")
