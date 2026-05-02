#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V119_Integration")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V119: 阈值特征证据整合")
message("========================================")

# --- 收集所有V114-V118结果 ---
files <- c(
  "03_results/V114_Bimodal/V114_bimodal_summary.csv",
  "03_results/V115_Cross_Model/V115_cross_model_bimodal.csv",
  "03_results/V116_Acute/V116_acute_bimodal.csv",
  "03_results/V117_Independent/V117_independent_bimodal.csv"
)

df_all <- NULL
for (f in files) {
  if (file.exists(f)) {
    d <- read.csv(f, stringsAsFactors=FALSE)
    df_all <- rbind(df_all, d)
    message("[LOAD] ", basename(f), ": ", nrow(d), " rows")
  } else {
    message("[SKIP] ", f, " not found")
  }
}

if (is.null(df_all) || nrow(df_all)==0) {
  message("[FAIL] No bimodal results found. Please run V114-V118 first.")
  quit(status=1)
}

# --- 统一叙事表 ---
df_all$threshold_signature <- ifelse(df_all$all_or_none > 70, "All-or-None (Threshold)", 
                                     ifelse(df_all$all_or_none > 50, "Partial Threshold", "Gradual"))
df_all$evidence_strength <- ifelse(df_all$is_bimodal==TRUE & df_all$all_or_none>70, "STRONG",
                                     ifelse(df_all$is_bimodal==TRUE, "MODERATE", "WEAK"))

write.csv(df_all, file.path(outdir,"V119_threshold_evidence_master_table.csv"), row.names=FALSE)
message("\n=== V119 阈值证据总表 ===")
print(df_all[,c("label","n","is_bimodal","zero_pct","all_or_none","threshold_signature","evidence_strength")])

# --- 可视化1: 全或无指数雷达 ---
p1 <- ggplot(df_all, aes(x=reorder(label, all_or_none), y=all_or_none, fill=evidence_strength)) +
  geom_bar(stat="identity") + coord_flip() +
  scale_fill_manual(values=c("STRONG"="#FDE725","MODERATE"="#35B779","WEAK"="#440154")) +
  geom_hline(yintercept=70, linetype="dashed", color="red") +
  annotate("text", x=1, y=75, label="Threshold Line (70%)", color="red", hjust=0) +
  labs(title="All-or-None Index Across All Models", x="", y="Zero% + Top10% (All-or-None Index)") +
  theme_minimal(base_size=10)

# --- 可视化2: 双峰性 vs 样本量 ---
p2 <- ggplot(df_all, aes(x=n, y=gap, color=evidence_strength, size=all_or_none)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values=c("STRONG"="#FDE725","MODERATE"="#35B779","WEAK"="#440154")) +
  labs(title="Bimodality Gap vs Sample Size", x="N cells/samples", y="Peak Gap") +
  theme_minimal()

# --- 可视化3: 零值比例比较 ---
p3 <- ggplot(df_all, aes(x=reorder(label, zero_pct), y=zero_pct, fill=threshold_signature)) +
  geom_bar(stat="identity") + coord_flip() +
  labs(title="NDUFB7 Silence Rate (Zero Expression)", x="", y="% Cells/Samples with NDUFB7=0") +
  theme_minimal()

# --- 组合图 ---
combined <- p1 + p2 + p3 + plot_layout(ncol=1, heights=c(2,1,1))
ggsave(file.path(outdir,"V119_threshold_evidence_figure.png"), combined, width=10, height=12, dpi=300)

# --- 生成Results段落框架 ---
cat("
=== Results Paragraph Framework (Threshold Signature) ===

Paragraph X: NDUFB7 Exhibits a Threshold-Characteristic Collapse Pattern

Across [N] independent datasets encompassing [list], NDUFB7 displayed a consistent
'all-or-none' expression pattern rather than gradual downregulation. In GSE183852
(n=220,752 cardiomyocyte nuclei), [X]% of cells showed complete NDUFB7 silence
(zero expression), while [Y]% maintained high expression, leaving only [Z]% in
an intermediate state — a distribution consistent with a threshold-driven
regulatory switch (bimodal mixture model, G=2, gap=[value]).

This pattern was [conserved/not conserved] across species and etiologies:
- Human DCM/ICM bulk: [all_or_none]%
- Mouse HFpEF: [all_or_none]%
- Embryonic development: [all_or_none]% ([gradual/threshold] as expected)
- Acute injury (DOX/STEMI): [all_or_none]%

The threshold signature suggests that NDUFB7 loss is not a passive consequence
of generalized metabolic decline but an active, switch-like cellular decision
point — potentially driven by [mechanism: e.g., epigenetic silencing,
microenvironmental hypoxia threshold, or proteostatic collapse].

", file=file.path(outdir,"V119_results_draft.txt"))

message("\n[DONE] V119 Integration complete")
message("  Master table: ", file.path(outdir,"V119_threshold_evidence_master_table.csv"))
message("  Figure: ", file.path(outdir,"V119_threshold_evidence_figure.png"))
message("  Results draft: ", file.path(outdir,"V119_results_draft.txt"))
message("\n[STRATEGY] 若跨模型中≥3个数据集 all_or_none>70% + is_bimodal=TRUE，则'阈值特征'可写入Results核心段落")
message("[STRATEGY] 若仅GSE183852满足，则降级为'在单细胞层面观察到阈值样分布'，不推广为跨模型普适规律")
