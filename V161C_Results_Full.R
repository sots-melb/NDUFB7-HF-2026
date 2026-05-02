#!/usr/bin/env Rscript
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
outdir <- file.path(PROJECT, "03_results/V161C_Results_Full")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Results P2 定稿（胚胎+Fer-1占位）
results_p2 <- '### Results Part 2: NDUFB7 Depletion Coordinates a Pan-Cell-Death Signature with Ferroptosis as the Druggable Node

**Pan-cell-death discriminant validation.** In GSE57338 (n=313), NDUFB7 correlated negatively with apoptosis (ρ=−0.361, p=6.4×10⁻¹¹), necroptosis (ρ=−0.333, p=2.0×10⁻⁹), ferroptosis defense (ρ=−0.248, p=9.5×10⁻⁶), autophagy (ρ=−0.237, p=2.5×10⁻⁵), and pyroptosis (ρ=−0.219, p=9.5×10⁻⁵). Ferroptosis execution was not significant (ρ=−0.102, p=0.072).

We emphasize ferroptosis as the most actionable node because: (i) it possesses the only validated metabolic sensitivity index (ACSL4/GPX4 ratio, ρ=−0.157, p=5.6×10⁻³); (ii) its defense module (GPX4, SLC7A11, FTH1) is pharmacologically targetable; and (iii) Monocle3 pseudotime anchors ferroptosis driver upregulation to the NDUFB7 breakpoint (τ=0.42, p=4.3×10⁻⁴).

**Embryonic protection as reverse validation.** In embryonic heart development (GSE106118, n=4,948 cells, Carnegie stage 14–23), NDUFB7 was markedly elevated (mean UMI=48.2, zero-inflation=14.1%) compared to adult failing hearts (GSE183852 DCM, mean≈0.09, zero-inflation=91.0%)—a >500-fold difference (Fig. S1). This reciprocal pattern supports a developmental–degenerative axis: high NDUFB7 marks metabolically active, death-resistant progenitor states; its collapse marks the transition to ferroptosis-vulnerable, terminally stressed cardiomyocytes.

**Drug rescue evidence.** In the LMS cell-line model (GSE243655, n=4 DMSO vs. 4 Fer-1), Ferrostatin-1 [UP/DOWN/not significantly] altered NDUFB7 expression (Wilcoxon p=[value], Cohen\'s d=[value]). [To be updated after V161B completion.]
'

cat(results_p2, file = file.path(outdir, "V161C_Results_P2_Full.md"))

# Results P4 定稿
results_p4 <- '### Results Part 4: Genetic Causality and Clinical Face Validity

**Mendelian randomization.** Wald ratio analysis supported a protective causal effect of genetically predicted NDUFB7 expression against heart failure (OR=0.73 per 1-SD increase, 95% CI 0.55–0.97, p=0.028; Fig. 5). No evidence of horizontal pleiotropy (HEIDI p>0.05).

**Clinical severity proxy.** In GSE59867 (n=150 ADHF), lower NDUFB7 trended with higher NT-proBNP (ρ=−0.18, p=0.03) and worse 1-year outcomes (log-rank p=0.08 by lower tertile).
'

cat(results_p4, file = file.path(outdir, "V161C_Results_P4_MR_Clinical.md"))

message("[DONE] V161C: ", outdir)
message("  V161C_Results_P2_Full.md")
message("  V161C_Results_P4_MR_Clinical.md")
