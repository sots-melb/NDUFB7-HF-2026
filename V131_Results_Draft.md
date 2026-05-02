# Results Paragraph 3-4 草稿（基于V113B/V120/V126/V127/V92/V84A/V128）

## Paragraph 3: Stepwise NDUFB7 Depletion Identifies a Metabolic Crisis Subpopulation

To characterize the dynamics of NDUFB7 loss in failing cardiomyocytes, we performed Monocle3 pseudotime trajectory analysis on GSE183852 single-nucleus data (n=503 CM). NDUFB7 expression followed a stepwise depletion pattern rather than gradual downregulation, with a sharp breakpoint at pseudotime 6.14 (early vs. mid mean: 4.37 vs. 3.41, p=4.3×10⁻⁴, permutation-validated p<0.001 across 1,000 randomizations). Notably, late-pseudotime cells exhibited partial compensatory rebound (mean=4.11), suggesting incomplete metabolic rescue in surviving cardiomyocytes.

This trajectory resolved into seven transcriptional subclusters, among which Cluster 3 (n=48 cells, 7.5% of CM) emerged as a distinct NDUFB7-silent population with 72.9% zero expression. Differential expression analysis confirmed that Cluster 3 was not a technical artifact: mitochondrial content was not elevated (excluding dead-cell contamination), and—critically—6 of 7 Complex I subunits were significantly downregulated (all p<0.05) while EndMT/fibrosis markers remained silent (0/7 upregulated), arguing against lineage transformation. The top downregulated genes included GAPDH and B2M, indicating a state of global metabolic suppression beyond isolated Complex I failure.

[Figure 3C: Monocle3 stepwise plot with breakpoint annotation]
[Figure 3D: Cluster 3 violin plot + DE heatmap]

---

## Paragraph 4: OXPHOS Collapse Drives Ferroptosis Defense Vulnerability

To test whether NDUFB7 loss propagates to multi-complex OXPHOS dysfunction, we compared Complex I–V expression between NDUFB7-low and NDUFB7-high CM. Complex I (log₂FC=−0.453, p=1.98×10⁻¹⁴), Complex III (−0.471, p=2.12×10⁻⁶), and Complex IV (−0.458, p=9.37×10⁻²⁶) showed coordinated collapse, whereas Complex II (nuclear-encoded, log₂FC=−0.066, p=0.697) was spared—consistent with its independent assembly and lower ROS sensitivity. This pattern supports a mitochondrial-specific crisis rather than generalized transcriptional shutdown.

Mechanistically, we quantified ferroptosis susceptibility using the ACSL4/GPX4 ratio index. NDUFB7-zero cells exhibited the highest ratio (25.96 vs. 24.15 in NDUFB7-high, Δ=+7.5%), indicating increased pro-ferroptotic drive relative to antioxidant defense. AUCell background-corrected scoring independently confirmed this direction (ferroptosis score: −0.421 in zero vs. −0.240 in high), whereas uncorrected sum-score was confounded by metabolic-active-cluster bias, underscoring the necessity of background normalization for low-expressed pathway genes.

Finally, CellChat analysis of CM–CM communication (n=3,742 significant pairs) identified ANXA1–FPR1 as the dominant axis (median probability=0.064 vs. background 3×10⁻⁴, p=6.7×10⁻⁵), with preferential signaling from NDUFB7-deficient toward NDUFB7-intermediate subclusters. ANXA1 is a classic pro-resolving mediator, suggesting that metabolically compromised CM actively emit anti-inflammatory "distress signals" to neighboring cells—a paradoxical rescue attempt by cells that are themselves unable to restore OXPHOS function.

[Figure 4A–D: Mechanism integration — stepwise→OXPHOS→ferroptosis→vicious cycle]
[Figure 5A: CellChat interaction heatmap]
[Figure 5D: ANXA1-FPR1 directionality]

---

## Honesty Notes for Methods/Discussion

1. **Monocle3 breakpoint**: Post-hoc segmentation; p-value not corrected for multiple window sizes. Permutation test (n=1,000) supports non-randomness but independent cohort validation is needed.
2. **Cluster 3 sample size**: n=48 is small; we acknowledge this as an exploratory subpopulation requiring validation in larger cohorts.
3. **CellChat**: In silico prediction tool; ANXA1–FPR1 directionality is computational inference, not experimental proof of intercellular communication.
4. **Ferroptosis**: "Susceptibility" not "execution" — direct lipid peroxidation evidence (4-HNE, C11-BODIPY) is absent and noted as a limitation.
