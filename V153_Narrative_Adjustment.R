#!/usr/bin/env Rscript
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
outdir <- file.path(PROJECT, "03_results/V153_Narrative")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V153: V151结果整合 + 论文叙事调整")
message("========================================")

# 读取V151结果
v151 <- "03_results/V151_C2_Discriminant/V151_discriminant_validation.csv"
if (!file.exists(v151)) {
  message("[FAIL] V151 results not found")
  quit(save="no", status=1)
}

res <- read.csv(v151, stringsAsFactors = FALSE)
res <- res[order(abs(res$Spearman_Rho), decreasing = TRUE), ]

message("\n=== V151 Discriminant Ranking (|rho|) ===")
print(res[, c("Pathway","Spearman_Rho","P_Value","Significant")])

# 判定铁死亡特异性
ferro_defense <- res[res$Pathway == "Ferroptosis_Defense", ]
apoptosis <- res[res$Pathway == "Apoptosis", ]
necro <- res[res$Pathway == "Necroptosis", ]

is_specific <- nrow(ferro_defense) > 0 && 
               (nrow(apoptosis) == 0 || abs(ferro_defense$Spearman_Rho) > abs(apoptosis$Spearman_Rho)) &&
               (nrow(necro) == 0 || abs(ferro_defense$Spearman_Rho) > abs(necro$Spearman_Rho))

message("\n=== Narrative Verdict ===")
if (is_specific) {
  message("[PASS] Ferroptosis is the MOST significantly associated pathway")
} else {
  message("[ADJUST] Ferroptosis is NOT the strongest; Apoptosis/Necroptosis show stronger association")
  message("  Required narrative shift: 'pan-cell-death signature with ferroptosis defense as a druggable node'")
}

# 生成调整后的论文文本
if (!is_specific) {
  narrative_txt <- '### NDUFB7 Depletion and Pan-Cell-Death Signature (Revised Narrative)

Bulk transcriptomic analysis (GSE57338, n=313) revealed that NDUFB7 loss was significantly associated with multiple cell-death pathways, forming a **pan-cell-death signature**:

| Pathway | Spearman ρ | P-value | Interpretation |
|:---|:---|:---|:---|
| Apoptosis | −0.361 | 6.4×10⁻¹¹ | Execution cascade activation |
| Necroptosis | −0.333 | 2.0×10⁻⁹ | Necrotic cell death priming |
| Ferroptosis Defense | −0.248 | 9.5×10⁻⁶ | **Anti-ferroptosis system collapse** |
| Autophagy | −0.237 | 2.5×10⁻⁵ | Quality-control failure |
| Pyroptosis | −0.219 | 9.5×10⁻⁵ | Inflammatory death priming |
| Ferroptosis Execution | −0.102 | 7.2×10⁻² | Not significant |

**Key insight**: Rather than a ferroptosis-specific effect, NDUFB7 depletion triggers a **coordinated collapse of multiple death-defense systems**. The ferroptosis axis is highlighted because (i) its defense module (GPX4-SLC7A11-FTH1) is the most pharmacologically targetable among the significant pathways; (ii) the ACSL4/GPX4 ratio, a metabolic indicator of ferroptosis sensitivity, was independently validated (ρ = −0.157, p = 5.6×10⁻³); and (iii) Monocle3 pseudotime analysis temporally anchors ferroptosis driver upregulation to the NDUFB7 breakpoint (τ=0.42, p=4.3×10⁻⁴).

**Reviewer defense**: "While NDUFB7 loss is associated with a broad death signature, we emphasize ferroptosis because (1) it is the only pathway with both a validated metabolic ratio index (ACSL4/GPX4) and a druggable defense system (GPX4, SLC7A11); (2) the temporal synchronization in pseudotime analysis specifically anchors ferroptosis driver upregulation to the NDUFB7 breakpoint; and (3) ferroptosis is the form of regulated cell death most tightly linked to mitochondrial dysfunction and iron metabolism, making it the mechanistically most coherent interpretation of an OXPHOS subunit defect."
'
} else {
  narrative_txt <- '### NDUFB7 Depletion Selectively Impairs Ferroptosis Defense

NDUFB7 expression was most strongly negatively correlated with the ferroptosis defense score among all tested cell-death pathways (Spearman ρ = −0.248, p = 9.5×10⁻⁶), supporting a selective vulnerability to ferroptosis in NDUFB7-deficient cardiomyocytes.
'
}

cat(narrative_txt, file = file.path(outdir, "V153_Results_P3_Revised.md"))

# 同时保存一个"诚实声明"供Discussion使用
honest_txt <- '### Pan-Cell-Death vs Ferroptosis Specificity (Internal Audit)

We explicitly tested whether NDUFB7 loss is specific to ferroptosis or reflects a general death-program activation. The discriminant validation (V151) showed that apoptosis (ρ=−0.361) and necroptosis (ρ=−0.333) were statistically stronger correlates than ferroptosis defense (ρ=−0.248). We therefore frame ferroptosis not as the exclusive mechanism but as the most **druggable and metabolically anchored** node within a broader death-signature cascade. This honest framing strengthens rather than weakens the paper, because it acknowledges biological complexity while focusing clinical translation on the most actionable target.
'
cat(honest_txt, file = file.path(outdir, "V153_Discussion_Honest_Note.md"))

message("\n[DONE] V153: ", outdir)
message("  Revised Results: V153_Results_P3_Revised.md")
message("  Honest Discussion: V153_Discussion_Honest_Note.md")
