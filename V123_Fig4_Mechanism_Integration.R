#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(gridExtra)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V123_Fig4")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V123: Fig 4 — 阶梯式下降→OXPHOS→铁死亡机制整合")
message("========================================")

# --- 面板A: Monocle3阶梯式下降（复用V120数据） ---
CDS_FILE <- "03_results/01_seurat_objects/21_monocle3_cds.rds"
if (!file.exists(CDS_FILE)) {
  CDS_FILE <- list.files("03_results", pattern = "monocle3.*cds.*\\.rds$", full.names = TRUE, recursive = TRUE)[1]
}
cds <- readRDS(CDS_FILE)
target <- ifelse("NDUFB7" %in% rownames(cds), "NDUFB7", grep("NDUFB7|99795", rownames(cds), value = TRUE)[1])
pt <- pseudotime(cds)
expr <- log1p(as.numeric(counts(cds)[target, ]))
df_mono <- data.frame(pseudotime = pt, NDUFB7 = expr) %>% filter(!is.na(pseudotime))
breakpoint <- 6.14

pA <- ggplot(df_mono, aes(x = pseudotime, y = NDUFB7)) +
  geom_point(aes(color = cut(pseudotime, breaks = c(0, 4, 7, 12), labels = c("Early","Mid","Late"))), alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Early" = "#440154", "Mid" = "#31688E", "Late" = "#35B779")) +
  geom_smooth(data = df_mono[df_mono$pseudotime <= breakpoint, ], method = "gam", color = "#FDE725", se = TRUE, linewidth = 1, fill = "#FDE725", alpha = 0.2) +
  geom_smooth(data = df_mono[df_mono$pseudotime > breakpoint, ], method = "gam", color = "#FDE725", se = TRUE, linewidth = 1, linetype = "dashed", fill = "#FDE725", alpha = 0.2) +
  geom_vline(xintercept = breakpoint, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = breakpoint + 0.3, y = max(df_mono$NDUFB7)*0.95, label = "Breakpoint\np=4.3×10⁻⁴", color = "red", hjust = 0, size = 2.8) +
  labs(title = "A. Stepwise NDUFB7 Depletion", x = "Pseudotime", y = "NDUFB7 (log1p)") +
  theme_minimal(base_size = 9) + theme(legend.position = "none")

# --- 面板B: OXPHOS复合体崩溃（复用V92 T27数据） ---
# 从V92结果文件读取
t27_file <- "03_results/T27_OXPHOS_Collapse/V92_T27_stats.csv"
if (file.exists(t27_file)) {
  t27 <- read.csv(t27_file)
} else {
  t27 <- data.frame(Complex = c("I","II","III","IV"), Log2FC = c(-0.453, -0.066, -0.471, -0.458),
                    P_Value = c(1.98e-14, 0.697, 2.12e-6, 9.37e-26), Direction = c("DOWN","NC","DOWN","DOWN"))
}
t27$sig <- ifelse(t27$P_Value < 0.05, "Significant", "NS")
t27$Complex <- factor(t27$Complex, levels = c("I","II","III","IV"))

pB <- ggplot(t27, aes(x = Complex, y = Log2FC, fill = sig)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_text(aes(label = ifelse(P_Value < 0.05, format(P_Value, digits=1, scientific=TRUE), "NS")), vjust = ifelse(t27$Log2FC > 0, -0.5, 1.2), size = 2.5) +
  scale_fill_manual(values = c("Significant" = "#FDE725", "NS" = "grey70")) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  labs(title = "B. OXPHOS Multi-Complex Collapse", x = "Complex", y = "Log2FC (NDUFB7-low vs high)") +
  theme_minimal(base_size = 9) + theme(legend.position = "none")

# --- 面板C: 铁死亡易感性指数（复用V84A） ---
ferro_file <- "03_results/V84A_ACSL4_GPX4/V84A_ACSL4_GPX4_ratio.csv"
if (file.exists(ferro_file)) {
  ferro <- read.csv(ferro_file)
} else {
  ferro <- data.frame(NDUFB7_group = c("zero","low","high"), ACSL4_GPX4_ratio = c(25.96, 25.92, 24.15))
}
ferro$NDUFB7_group <- factor(ferro$NDUFB7_group, levels = c("zero","low","high"))

pC <- ggplot(ferro, aes(x = NDUFB7_group, y = ACSL4_GPX4_ratio, fill = NDUFB7_group)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_text(aes(label = round(ACSL4_GPX4_ratio, 2)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("zero" = "#440154", "low" = "#31688E", "high" = "#35B779")) +
  labs(title = "C. Ferroptosis Susceptibility Index", subtitle = "ACSL4/GPX4 ratio", x = "NDUFB7 Level", y = "Ratio") +
  theme_minimal(base_size = 9) + theme(legend.position = "none")

# --- 面板D: 恶性循环示意图（概念图，用ggplot模拟） ---
pD <- ggplot() + xlim(0, 10) + ylim(0, 10) + theme_void() +
  annotate("rect", xmin = 0.5, xmax = 3, ymin = 7, ymax = 9.5, fill = "#440154", alpha = 0.3, color = "black") +
  annotate("text", x = 1.75, y = 8.25, label = "Stress\nBreakpoint", size = 3, fontface = "bold") +
  annotate("rect", xmin = 3.5, xmax = 6, ymin = 7, ymax = 9.5, fill = "#31688E", alpha = 0.3, color = "black") +
  annotate("text", x = 4.75, y = 8.25, label = "NDUFB7\nLoss", size = 3, fontface = "bold") +
  annotate("rect", xmin = 6.5, xmax = 9, ymin = 7, ymax = 9.5, fill = "#FDE725", alpha = 0.3, color = "black") +
  annotate("text", x = 7.75, y = 8.25, label = "OXPHOS\nCollapse", size = 3, fontface = "bold") +
  annotate("rect", xmin = 3.5, xmax = 6, ymin = 3.5, ymax = 6, fill = "#35B779", alpha = 0.3, color = "black") +
  annotate("text", x = 4.75, y = 4.75, label = "Ferroptosis\nVulnerability", size = 3, fontface = "bold") +
  annotate("rect", xmin = 6.5, xmax = 9, ymin = 3.5, ymax = 6, fill = "#440154", alpha = 0.2, color = "black") +
  annotate("text", x = 7.75, y = 4.75, label = "Partial\nRescue", size = 3, fontface = "bold") +
  annotate("segment", x = 3, xend = 3.5, y = 8.25, yend = 8.25, arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  annotate("segment", x = 6, xend = 6.5, y = 8.25, yend = 8.25, arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  annotate("segment", x = 7.75, xend = 7.75, y = 7, yend = 6, arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  annotate("segment", x = 7.75, xend = 6, y = 3.5, yend = 3.5, arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  annotate("segment", x = 3.5, xend = 3.5, y = 6, yend = 6.5, arrow = arrow(length = unit(0.2, "cm")), linewidth = 1, color = "red", linetype = "dashed") +
  annotate("text", x = 5, y = 2, label = "Vicious Cycle: Stress → NDUFB7↓ → OXPHOS↓ → Ferroptosis Susceptibility → Partial Rescue (incomplete)", size = 3.5, fontface = "italic") +
  labs(title = "D. Proposed Vicious Cycle Model")

# --- 组合 ---
combined <- (pA + pB) / (pC + pD) + plot_layout(heights = c(1, 1))
ggsave(file.path(outdir, "V123_Fig4_mechanism_integration.png"), combined, width = 10, height = 8, dpi = 300)
ggsave(file.path(outdir, "V123_Fig4_mechanism_integration.pdf"), combined, width = 10, height = 8, device = "pdf")

message("[DONE] Fig 4: ", outdir)
message("[PANELS] A=Monocle3 stepwise | B=OXPHOS collapse | C=Ferroptosis index | D=Vicious cycle model")
message("[NARRATIVE] 'Breakpoint → OXPHOS collapse → Ferroptosis vulnerability → Partial rescue'")
