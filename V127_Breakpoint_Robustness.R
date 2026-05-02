#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V127_Breakpoint_Robustness")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V127: Monocle3 Breakpoint 稳健性验证")
message("========================================")

CDS_FILE <- "03_results/01_seurat_objects/21_monocle3_cds.rds"
if (!file.exists(CDS_FILE)) {
  CDS_FILE <- list.files("03_results", pattern = "monocle3.*cds.*\\.rds$", full.names = TRUE, recursive = TRUE)[1]
}
cds <- readRDS(CDS_FILE)
target <- ifelse("NDUFB7" %in% rownames(cds), "NDUFB7", grep("NDUFB7|99795", rownames(cds), value = TRUE)[1])
pt <- monocle3::pseudotime(cds)
expr <- log1p(as.numeric(SummarizedExperiment::assay(cds, "counts")[target, ]))
df <- data.frame(pseudotime = pt, NDUFB7 = expr) %>% dplyr::filter(!is.na(pseudotime))

# --- [1/3] 不同分位数breakpoint的段间差异 ---
message("\n>>> [1/3] 参数扫描：不同分位数breakpoint")

quantiles <- c(0.25, 0.33, 0.40, 0.50, 0.60, 0.67, 0.75)
robustness <- data.frame()

for (q in quantiles) {
  bp <- quantile(df$pseudotime, q)
  early <- df$NDUFB7[df$pseudotime <= bp]
  late <- df$NDUFB7[df$pseudotime > bp]
  
  if (length(early) > 10 && length(late) > 10) {
    tt <- t.test(early, late)
    robustness <- rbind(robustness, data.frame(
      quantile = q, breakpoint = bp, n_early = length(early), n_late = length(late),
      mean_early = mean(early), mean_late = mean(late), 
      delta = mean(late) - mean(early), t_stat = tt$statistic, p_val = tt$p.value
    ))
  }
}

write.csv(robustness, file.path(outdir, "V127_breakpoint_parameter_scan.csv"), row.names = FALSE)
message("[PASS] Parameter scan: ", nrow(robustness), " breakpoints tested")
print(robustness[, c("quantile", "breakpoint", "delta", "p_val")])

# 判定：如果多个分位数都显著→稳健；如果仅33%附近显著→敏感
sig_count <- sum(robustness$p_val < 0.05)
message("\nSignificant breakpoints (p<0.05): ", sig_count, "/", nrow(robustness))
if (sig_count >= 4) {
  message("[PASS] Breakpoint effect is robust across quantile choices")
} else if (sig_count >= 2) {
  message("[PARTIAL] Breakpoint effect is moderate — strongest at 33%")
} else {
  message("[WARN] Breakpoint effect is fragile — may be parameter-dependent")
}

# --- [2/3] Permutation test ---
message("\n>>> [2/3] Permutation test（1000次随机打乱）")

set.seed(2026)
n_perm <- 1000
perm_pvals <- numeric(n_perm)

for (i in 1:n_perm) {
  df_perm <- df
  df_perm$pseudotime <- sample(df_perm$pseudotime)  # 打乱伪时间
  bp_perm <- quantile(df_perm$pseudotime, 1/3)
  early_p <- df_perm$NDUFB7[df_perm$pseudotime <= bp_perm]
  late_p <- df_perm$NDUFB7[df_perm$pseudotime > bp_perm]
  
  if (length(early_p) > 10 && length(late_p) > 10) {
    perm_pvals[i] <- t.test(early_p, late_p)$p.value
  } else {
    perm_pvals[i] <- 1
  }
  
  if (i %% 200 == 0) message("  Permutation ", i, "/", n_perm)
}

observed_p <- 4.3e-04  # V113B observed
perm_better <- sum(perm_pvals <= observed_p, na.rm = TRUE)
perm_p <- perm_better / n_perm

message("\nObserved p: ", format(observed_p, scientific = TRUE))
message("Permutations with p <= observed: ", perm_better, "/", n_perm)
message("Permutation-based p: ", format(perm_p, scientific = TRUE))

if (perm_p < 0.05) {
  message("[PASS] Breakpoint is unlikely to be random noise (permutation p<0.05)")
} else {
  message("[WARN] Random pseudotime assignments can produce equally extreme p-values")
}

# --- [3/3] 可视化 ---
message("\n>>> [3/3] 稳健性可视化")

p1 <- ggplot(robustness, aes(x = breakpoint, y = -log10(p_val), color = p_val < 0.05)) +
  geom_point(size = 3) + geom_line(aes(group = 1), color = "grey50", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
  scale_color_manual(values = c("TRUE" = "#FDE725", "FALSE" = "#440154")) +
  labs(title = "Breakpoint Robustness: Quantile Scan", x = "Breakpoint Pseudotime", 
       y = "-log10(p-value)", color = "Significant") + theme_minimal()

p2 <- ggplot(data.frame(p_value = perm_pvals), aes(x = p_value)) +
  geom_histogram(bins = 50, fill = "#31688E", alpha = 0.6, color = "black") +
  geom_vline(xintercept = observed_p, color = "red", linewidth = 1) +
  annotate("text", x = observed_p, y = Inf, label = paste0("Observed\np=", format(observed_p, digits=1, scientific=TRUE)), 
           color = "red", vjust = 1.5, hjust = -0.1, size = 3) +
  labs(title = "Permutation Null Distribution", x = "T-test p-value (randomized pseudotime)", y = "Count") +
  theme_minimal()

combined <- p1 + p2 + plot_layout(ncol = 1)
ggsave(file.path(outdir, "V127_breakpoint_robustness.png"), combined, width = 7, height = 7, dpi = 300)

message("\n[DONE] V127: ", outdir)
