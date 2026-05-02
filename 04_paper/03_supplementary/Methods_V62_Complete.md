# Methods (V62 Complete, ~1200 words)

## Study Design

Integrative analysis combining spatial transcriptomics (Visium V1), multi-platform bulk RNA validation (5 cohorts, n=431), weighted gene co-expression network analysis (WGCNA), single-cell reference analysis (2 datasets, n=230,746 cells), protein-level validation, eQTL analysis, and in silico knockout pathway simulation.

## Spatial Transcriptomics

Kuppe Visium V1 (Zenodo 6578047, 6 hearts, 19,525 spots). Zones: FZ (n=2), IZ (n=3), IZ_BZ (n=1). NDUFB7 positivity compared by Fisher's exact test. Linear mixed-effects model (LMM) with sample ID as random effect applied to address spatial pseudoreplication (confirmed by Moran's I analysis). Published cell2location deconvolution used for cell-type composition inference.

## Bulk RNA Validation

Five cohorts: (1) GSE57338 Affy 1.1 ST (n=313); (2) GSE116250 RNA-seq RPKM (n=64, **length-bias artifact flagged**); (3) GSE55296 RNA-seq counts (n=37); (4) GDS4772 Affy U133A (n=17); (5) PXD010154 proteomics (19 fractions). Effect sizes: Cohen's d with 95% CI. Forest plot for cross-platform consistency.

## Single-Cell Analysis

GSE109816 (normal left atrium, n=9,994 cells): NCM NDUFB7 median=4 vs CM median=1 (p=4.4×10⁻⁸²). GSE183852 (donor/DCM left ventricle, n=220,752 nuclei): Donor CM mean=0.205 vs DCM CM mean=0.140 (p=5×10⁻²⁵); Donor FB mean=0.134 vs DCM FB mean=0.085 (p=1.4×10⁻⁵³). Cell types from published annotations (Names column). **Zero-inflation declared**: 87-92% zero values across all cell types.

## eQTL Analysis

eQTLGen Phase I (blood, n=31,684): rs11085898 (T allele increases NDUFB7, p=2.6×10⁻⁹). GTEx v11 Heart LV (n=387): rs8103021 (T allele decreases NDUFB7, p=1.89×10⁻⁶). **Opposite directional effects** highlight tissue-specific genetic regulation.

## In Silico Knockout

PROGENy pathway simulation: NDUFB7 bottom 25% vs top 25% (GSE57338, n=158). OXPHOS impairment (p=1.13×10⁻¹³), ROS defense attenuation (p=6.06×10⁻⁴), fibrosis activation (p=3.94×10⁻⁶). **Explicitly labeled as predictive, not causal**.

## Statistical Analysis

R v4.3.1. Two-sided tests. P < 0.05 significant. Benjamini-Hochberg FDR correction applied. Effect sizes reported alongside p-values. Code: https://github.com/sots-melb/NDUFB7-HF-2026 (MIT License).
