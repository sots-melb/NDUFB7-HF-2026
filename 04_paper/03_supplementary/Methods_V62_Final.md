# Methods (V62 Final, ~1500 words)

## Study Design and Data Sources

This integrative study combines spatial transcriptomics (Visium V1), multi-platform bulk RNA validation (5 cohorts, n=431), weighted gene co-expression network analysis (WGCNA), single-cell/nucleus reference analysis (2 datasets, n=230,746 cells), protein-level validation, expression quantitative trait locus (eQTL) analysis, and in silico knockout pathway simulation. All data were obtained from publicly available repositories. The study was designed to characterize the spatial dynamics, etiology-specific regulation, and cellular origin of NDUFB7 (NADH:ubiquinone oxidoreductase subunit B7, 137 amino acids), an accessory subunit of mitochondrial Complex I.

## Spatial Transcriptomics Analysis

Spatial transcriptomics data were obtained from Kuppe et al. (Zenodo 6578047, Nature 2022), comprising six post-myocardial infarction human heart slices (19,525 spots) profiled by 10x Genomics Visium V1. Tissue zones were annotated as fibrotic zone (FZ, n=2 hearts), ischemic zone (IZ, n=3), and IZ-border zone (IZ_BZ, n=1). NDUFB7 spot-level expression was extracted from normalized UMI counts. Zone-level positivity rates (spots with NDUFB7 > 0) were compared by Fisher's exact test.

To address spatial pseudoreplication, we applied **Moran's I test** (spdep R package) to quantify spatial autocorrelation in NDUFB7 expression. Significant positive autocorrelation (I > 0, p < 0.001) confirmed non-independence of spots within each slice. We therefore applied a **linear mixed-effects model (LMM)** with sample ID as random effect (lme4 R package, v1.1-35.1), which provides conservative estimates appropriate for spatially dependent data. Model comparison (AIC/BIC) confirmed LMM superiority over fixed-effects-only models (ΔAIC > 2).

Published cell2location deconvolution (Kuppe et al., Nature 2022) was leveraged to infer cell-type composition of Visium spots. NDUFB7-high (top quartile) and NDUFB7-low (bottom quartile) spots were compared for cardiomyocyte (CM) and fibroblast (FB) enrichment.

## Multi-Platform Bulk RNA Validation

Five independent cohorts were analyzed: (1) GSE57338 (Affymetrix Human Gene 1.1 ST Array, n=313, HF vs non-failing); (2) GSE116250 (Illumina HiSeq 2000, RNA-seq RPKM, n=64, DCM vs non-failing); (3) GSE55296 (Illumina HiSeq 2000, RNA-seq counts, n=37, DCM vs ischemic cardiomyopathy); (4) GDS4772 (Affymetrix U133A 2.0, n=17, DCM vs normal); (5) PXD010154 (LC-MS/MS proteomics, 19 fractions). NDUFB7 expression was extracted by platform-specific probe mapping (Affymetrix) or Ensembl ID alignment (RNA-seq). Effect sizes were computed as Cohen's d with 95% confidence intervals. Cross-platform consistency was visualized by forest plot.

**Methodological caution**: GSE116250 RPKM values may be inflated by length-normalization bias (NDUFB7 CDS = 411 bp, 137 aa), as RPKM divides by transcript length. This artifact was flagged and discussed; conclusions prioritize count-based platforms (GSE55296).

## Single-Cell and Single-Nucleus Analysis

Two datasets were analyzed: (1) GSE109816 (normal human left atrium, n=9,994 cells, Litviňuková et al., Nature 2020); (2) GSE183852 (donor and DCM human left ventricle, n=220,752 nuclei, Koenig et al., Nature 2022). Cell-type annotations were obtained from published metadata (Names column: Cardiomyocytes, Fibroblasts, Endothelium, etc.). NDUFB7 UMI counts were extracted per cell type and condition (Donor vs DCM). **Wilcoxon rank-sum test** was used for group comparisons.

**Zero-inflation declaration**: NDUFB7 showed extreme sparse expression (87-92% zero values across all cell types), with median = 0 in most populations. This precludes conventional parametric differential expression analysis. We report zero-inflation rates, medians, and means alongside Wilcoxon p-values to transparently characterize expression distribution.

## Expression Quantitative Trait Locus (eQTL) Analysis

Cis-eQTL data were obtained from eQTLGen Phase I (blood, n = 31,684) and GTEx v11 (heart left ventricle, n = 387). Lead SNPs within ±1 Mb of NDUFB7 were extracted. **Critical finding**: eQTLGen (rs11085898, T allele increases NDUFB7, p = 2.6×10⁻⁹) and GTEx Heart-LV (rs8103021, T allele decreases NDUFB7, p = 1.89×10⁻⁶) showed opposite directional effects, highlighting tissue-specific genetic regulation. Ensembl ID verification confirmed GTEx v11 uses ENSG00000099795.7 (not ENSG00000167996).

## In Silico Knockout and Pathway Simulation

PROGENy (R package, v1.20.0) was used to infer pathway activity from gene expression. GSE57338 samples were stratified by NDUFB7 expression (bottom 25% as "virtual knockout" vs top 25% as "virtual wild-type"). Pathway scores (OXPHOS, ROS, Hypoxia, TGF-β, etc.) were compared by Wilcoxon test. **Explicit caveat**: This is an in silico correlation-based prediction, not experimental knockout. Results are presented as mechanistic hypotheses requiring functional validation.

## Statistical Analysis

All statistical tests were two-sided. P < 0.05 was considered significant. **Benjamini-Hochberg false discovery rate (FDR) correction** was applied across all 15 reported tests (Supplementary Table S1). Effect sizes (Cohen's d) were reported alongside p-values. Analyses were performed in R v4.3.1. Code is available at https://github.com/sots-melb/NDUFB7-HF-2026 (MIT License).
