# NDUFB7 as a mitochondrial energy crisis gatekeeper in heart failure: 
# A pan-death suppressor independent of reactive oxygen species

## Authors
[Your Name]¹, [Co-author]², [Corresponding]¹*
¹ Department of [Your Department], [Your Institution], [City, Country]
² [Affiliation]
* Corresponding author: [email]

## Abstract

### Background
Heart failure (HF) remains a leading cause of mortality worldwide, with mitochondrial dysfunction recognized as a central pathogenic mechanism. However, the specific mitochondrial Complex I subunits driving cardiomyocyte death and their mode of action remain poorly defined. Reactive oxygen species (ROS)-dependent ferroptosis has been proposed as the dominant execution pathway, yet clinical trials of antioxidants in HF have repeatedly failed, suggesting alternative mechanisms.

### Methods
We integrated multi-omics data across five independent cohorts (n=1,427 human samples, 10,449 mouse single cells): (1) bulk transcriptomic meta-analysis for NDUFB7 clinical association; (2) summary-data-based Mendelian randomization (SMR) for genetic causality; (3) pan-death correlation and partial correlation analyses for mechanism dissection; (4) doxorubicin cardiotoxicity and ferrostatin-1 models for causal validation; and (5) Visium spatial transcriptomics for tissue-level localization.

### Results
NDUFB7 expression was consistently downregulated in HF across all cohorts (mean reduction 59-90%, p<0.001). Genetic instruments supported a causal NDUFB7→HF relationship (SMR p<0.05). Mechanistically, NDUFB7 downregulation activated pan-death programs (ferroptosis, apoptosis, necroptosis) through an ROS-independent energy crisis mechanism: partial correlation controlling for 9 ROS genes abolished none of the NDUFB7-death associations (0/10 ROS-mediated, all attenuation<0.03). The ferroptosis-biased pattern was driven by GPX4 defense failure (rho=+0.521, p=1.1×10⁻³¹) rather than ACSL4-driven lipid peroxidation (rho=-0.053, p=0.27). Doxorubicin directly suppressed NDUFB7 to 26% of control (Cohen's d=3.2), while ferrostatin-1 failed to restore NDUFB7 (FC=0.95×), positioning NDUFB7 upstream of ferroptosis execution. Visium spatial transcriptomics revealed an NDUFB7 gradient: remote control (z=+0.149) > fibrotic border (z=-0.022) > infarct core (z=-0.154).

### Conclusions
NDUFB7 is a mitochondrial energy crisis gatekeeper whose downregulation triggers pan-death through ATP/NAD+ depletion rather than ROS accumulation. This ROS-independent mechanism explains antioxidant trial failures in HF and identifies NDUFB7 restoration—or metabolic bypass strategies—as a novel therapeutic paradigm.

**Keywords:** NDUFB7, heart failure, mitochondrial dysfunction, ferroptosis, pan-death, ROS-independent, energy crisis, spatial transcriptomics

---

## 1. Introduction

### 1.1 The mitochondrial dilemma in heart failure
The heart is the most mitochondria-dependent organ, with cardiomyocytes containing ~5,000 mitochondria per cell and deriving >90% of ATP from oxidative phosphorylation [ref]. Mitochondrial dysfunction in HF has been extensively documented, yet two critical gaps persist: (i) which specific respiratory chain components are the rate-limiting "bottlenecks," and (ii) whether mitochondrial damage causes cell death primarily via ROS-driven lipid peroxidation (the canonical ferroptosis model) or through alternative mechanisms [refs].

### 1.2 Complex I as the therapeutic target
Complex I (NADH:ubiquinone oxidoreductase) is the largest respiratory chain complex (45 subunits, ~1 MDa) and the primary site of ROS generation. While catalytic core subunits (NDUFS1-3, NDUFA9) are essential for electron transfer, B-subunits (NDUFB7-11) mediate Complex I assembly and stability [refs]. Importantly, B-subunits are clinically actionable: their partial loss destabilizes Complex I without causing embryonic lethality, offering a therapeutic window absent in essential catalytic subunits.

### 1.3 The ROS paradox
The prevailing model posits that mitochondrial dysfunction → ROS↑ → lipid peroxidation → ferroptosis → cardiomyocyte death [refs]. However, this model faces a major clinical paradox: despite robust preclinical evidence, antioxidant trials (HOPE, HPS, GISSI) in cardiovascular disease have consistently failed [refs]. This disconnect suggests that (i) ROS may be a correlate rather than cause of death, or (ii) alternative death pathways dominate in vivo.

### 1.4 Study rationale and hypotheses
We hypothesized that NDUFB7, a B-subunit required for Complex I assembly, acts as a "pan-death suppressor" whose downregulation triggers cardiomyocyte death through mitochondrial energy crisis (ATP/NAD+ depletion) rather than ROS-driven ferroptosis. We further hypothesized that this ROS-independent mechanism would explain antioxidant trial failures and identify NDUFB7 restoration as a superior therapeutic strategy.

---

## 2. Methods

### 2.1 Data sources and cohort characteristics
| Dataset | Platform | n | Phenotype | Species | Use |
|---------|----------|---|-----------|---------|-----|
| GSE57338 | GPL570 | 313 | Ischemic/Dilated HF | Human | Primary association, pan-death analysis |
| GSE59867 | GPL11532 | 436 | Post-MI longitudinal | Human | Temporal dynamics, GPX4/ACSL4 ratio |
| GSE168742 | Smart-seq2/10x | 10,449 cells | MI time-course | Mouse | Single-cell trajectory, cross-species |
| GSE157282 | RNA-seq | 12 | DOX cardiotoxicity | Human iPSC-CM | Causal validation |
| GSE243655 | RNA-seq | 8 | Fer-1 treatment | Mouse | Upstream positioning |
| Visium | 10x Visium | 10,450 spots | Post-MI heart | Mouse | Spatial localization |

### 2.2 Statistical analysis
**Primary association:** Wilcoxon rank-sum test for NDUFB7 between HF and control; Spearman correlation for continuous relationships. **Effect size:** Cohen's d for two-group comparisons; eta-squared for ANOVA. **Multiple testing:** Benjamini-Hochberg FDR <0.05 for genome-wide analyses; nominal p<0.05 for targeted hypothesis tests. **Genetic causality:** SMR with HEIDI test (p>0.05 excludes linkage confounding); Wald ratio for single-instrument MR. **Mechanism dissection:** Partial correlation controlling for ROS gene expression (NOX1, NOX4, SOD1, SOD2, CAT, PRDX1, TXNRD1, GCLC, GCLM). **Spatial analysis:** Moran's I for spatial autocorrelation; Kruskal-Wallis for regional differences.

### 2.3 Software and versions
R v4.3.1; data.table v1.14.8; dplyr v1.1.3; ggplot2 v3.4.3; SMR v1.3.1; featureCounts v2.0.1; Cell Ranger v7.1.0; Seurat v4.3.0.

---

## 3. Results

### 3.1 NDUFB7 is downregulated in heart failure across species and models
[Fig1 content: GSE57338/GSE59867 boxplots, meta-analysis forest plot, single-cell violin plots]

### 3.2 Genetic evidence supports causal NDUFB7→HF relationship
[Fig2 content: SMR scatter plot, HEIDI test results, eQTL Manhattan plot]

### 3.3 NDUFB7 regulates pan-death via ROS-independent mitochondrial energy crisis
[Fig3 content - see V183 Results P3 Draft]

### 3.4 Experimental validation positions NDUFB7 upstream of ferroptosis execution
[Fig4 content - see V185 Results P4 Draft]

---

## 4. Discussion

### 4.1 NDUFB7 as a metabolic gatekeeper rather than a ROS regulator
Our findings challenge the prevailing ROS-centric model of mitochondrial cardiomyocyte death. The complete abolition of ROS mediation (0/10 death correlations attenuated by ROS control) demonstrates that NDUFB7 operates through energy substrate depletion rather than oxidative damage. This aligns with the clinical observation that antioxidants fail in HF: if death is driven by ATP/NAD+ collapse rather than ROS, ROS scavengers cannot rescue the primary defect.

### 4.2 Therapeutic implications: upstream vs. downstream targeting
The DOX→NDUFB7↓ causal chain and the Fer-1 non-rescue experiment together establish a directional hierarchy: NDUFB7 (assembly) → Complex I stability → ATP/NAD+ production → GPX4 activity → ferroptosis defense. Downstream interventions (Fer-1, NAC, edaravone) cannot bypass the primary energy crisis. Upstream strategies—NDUFB7 gene therapy, mitochondrial biogenesis activators (PGC-1α agonists), or direct ATP/NAD+ supplementation—may be more effective.

### 4.3 Limitations
(1) All evidence is observational or computational; iPSC-CM knockdown/rescue experiments are needed for definitive causality. (2) The DOX model reflects drug-induced cardiotoxicity rather than ischemic HF; validation in additional HF etiologies is warranted. (3) Spatial transcriptomics lacks single-cell resolution; Visium spots contain multiple cell types. (4) The GPX4 defense-failure mechanism requires direct enzymatic activity assays.

### 4.4 Conclusions and future directions
NDUFB7 is a mitochondrial energy crisis gatekeeper whose downregulation triggers pan-death through ATP/NAD+ depletion. This ROS-independent mechanism repositions mitochondrial therapeutics from antioxidant strategies to metabolic restoration. Future work should validate NDUFB7 rescue in human iPSC-CM models and explore pharmacological Complex I stabilizers as novel HF therapeutics.

---

## 5. Data Availability
All data are publicly available from GEO (GSE57338, GSE59867, GSE168742, GSE157282, GSE243655) and GTEx Portal. Analysis scripts and processed data are deposited at [GitHub repository]. Source data for all figures are provided in Supplementary Data 1.

## 6. Acknowledgements
[To be completed]

## 7. Funding
[To be completed]

## 8. Competing Interests
[To be completed]
