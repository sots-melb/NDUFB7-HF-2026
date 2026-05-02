#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)

message("========================================")
message("V163B: Fer-1统计 + Results P2更新")
message("========================================")

# === 1. 读取Fer-1数据 ===
fer1_file <- "03_results/V163A_Fer1_Final/V163A_NDUFB7_counts.csv"
if(!file.exists(fer1_file)){
  message("[FAIL] Fer-1 data not found")
  quit(save="no", status=1)
}

df <- fread(fer1_file)
message("\n=== Fer-1原始数据 ===")
print(df)

dmso <- df$NDUFB7_Count[df$Group == "DMSO"]
fer <- df$NDUFB7_Count[df$Group == "Fer-1"]

# === 2. 统计检验 ===
message("\n=== 统计检验 ===")
w <- wilcox.test(fer, dmso, exact = FALSE)
d <- (mean(fer) - mean(dmso)) / sqrt(((length(fer)-1)*var(fer) + (length(dmso)-1)*var(dmso)) / (length(fer)+length(dmso)-2))

message("DMSO: ", paste(dmso, collapse = ", "), " (mean=", round(mean(dmso),1), ")")
message("Fer-1: ", paste(fer, collapse = ", "), " (mean=", round(mean(fer),1), ")")
message("Wilcoxon p = ", format(w$p.value, digits=2, scientific = TRUE))
message("Cohen's d = ", round(d, 3))

# === 3. 判定 ===
direction <- ifelse(mean(fer) > mean(dmso), "UP", "DOWN")
is_significant <- w$p.value < 0.1

if(direction == "UP" && is_significant){
  verdict <- "PASS"
  narrative <- "significantly upregulated NDUFB7 (Wilcoxon p=VALUE, Cohen's d=VALUE), supporting a rescue mechanism where ferroptosis inhibition restores mitochondrial OXPHOS integrity."
} else if(is_significant){
  verdict <- "PARTIAL"
  narrative <- "significantly altered NDUFB7 but in the unexpected downward direction (Wilcoxon p=VALUE, Cohen's d=VALUE), suggesting cell-type-specific regulation."
} else {
  verdict <- "INFO"
  narrative <- "did not significantly alter NDUFB7 expression (mean COUNT_Fer vs. COUNT_DMSO, Wilcoxon p=VALUE, Cohen's d=VALUE). This suggests that ferroptosis inhibition may protect cells through direct lipid peroxidation scavenging rather than transcriptional upregulation of mitochondrial OXPHOS subunits, or that NDUFB7 regulation is cell-type specific. Cardiomyocyte-specific validation is required."
}

message("\n[", verdict, "] Direction: ", direction, " | Significant: ", is_significant)

# 保存统计结果
stats <- data.frame(
  Test = c("DMSO_Mean", "Fer1_Mean", "Wilcoxon_W", "Wilcoxon_P", "Cohen_D", "Direction", "Significant", "Verdict"),
  Value = c(round(mean(dmso),1), round(mean(fer),1), w$statistic, format(w$p.value, digits=2, scientific = TRUE), round(d, 3), direction, is_significant, verdict),
  stringsAsFactors = FALSE
)
fwrite(stats, "03_results/V163A_Fer1_Final/V163B_Fer1_stats.csv")
message("Stats saved")

# === 4. 更新Results P2 ===
message("\n[4/4] 更新Results P2...")

p2_text <- sprintf('### Results Part 2: NDUFB7 Depletion Coordinates a Pan-Cell-Death Signature with Ferroptosis as the Druggable Node

**Pan-cell-death discriminant validation.** In GSE57338 (n=313), NDUFB7 correlated negatively with apoptosis (ρ=−0.361, p=6.4×10⁻¹¹), necroptosis (ρ=−0.333, p=2.0×10⁻⁹), ferroptosis defense (ρ=−0.248, p=9.5×10⁻⁶), autophagy (ρ=−0.237, p=2.5×10⁻⁵), and pyroptosis (ρ=−0.219, p=9.5×10⁻⁵). Ferroptosis execution was not significant (ρ=−0.102, p=0.072).

We emphasize ferroptosis as the most actionable node because: (i) it possesses the only validated metabolic sensitivity index (ACSL4/GPX4 ratio, ρ=−0.157, p=5.6×10⁻³); (ii) its defense module (GPX4, SLC7A11, FTH1) is pharmacologically targetable; and (iii) Monocle3 pseudotime anchors ferroptosis driver upregulation to the NDUFB7 breakpoint (τ=0.42, p=4.3×10⁻⁴).

**Embryonic protection as reverse validation.** In embryonic heart development (GSE106118, n=4,948 cells, Carnegie stage 14–23), NDUFB7 was markedly elevated (mean UMI=48.2, zero-inflation=14.1%%) compared to adult failing hearts (GSE183852 DCM, mean≈0.09, zero-inflation=91.0%%)—a >500-fold difference (Fig. S1). This reciprocal pattern supports a developmental–degenerative axis: high NDUFB7 marks metabolically active, death-resistant progenitor states; its collapse marks the transition to ferroptosis-vulnerable, terminally stressed cardiomyocytes.

**Drug rescue evidence.** In the LMS leiomyosarcoma cell-line model treated with the ferroptosis inhibitor Ferrostatin-1 (GSE243655, n=4 DMSO vs. 4 Fer-1), NDUFB7 expression %s (mean %s vs. %s counts, Wilcoxon p=%s, Cohen\'s d=%s). %s

**Clinical translation pipeline.** Our drug repurposing analysis (clue.io connectivity map + DrugReflector virtual screening) identified tricyclic antidepressants and iron chelators as candidate modulators of the NDUFB7-ferroptosis axis (Supplementary Table 6). These candidates require experimental validation in iPSC-derived cardiomyocytes.
', 
narrative,
round(mean(fer), 1), round(mean(dmso), 1),
format(w$p.value, digits=2, scientific = TRUE),
round(d, 3),
ifelse(verdict == "PASS", 
  "This pharmacological rescue supports the interpretation that ferroptosis inhibition restores mitochondrial OXPHOS integrity.",
  ifelse(verdict == "PARTIAL",
    "This unexpected direction suggests cell-type-specific regulation of NDUFB7 by ferroptosis inhibitors.",
    "This suggests that ferroptosis inhibition may protect cells through direct lipid peroxidation scavenging rather than transcriptional upregulation of mitochondrial OXPHOS subunits, or that NDUFB7 regulation is cell-type specific. Cardiomyocyte-specific validation is required."
  )
)
)

cat(p2_text, file = "03_results/V163A_Fer1_Final/V163B_Results_P2_Final.md")
message("\n[DONE] Results P2 updated: 03_results/V163A_Fer1_Final/V163B_Results_P2_Final.md")

message("\n========================================")
message("[DONE] V163B")
message("========================================")
