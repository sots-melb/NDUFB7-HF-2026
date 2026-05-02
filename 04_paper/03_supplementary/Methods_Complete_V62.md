# Methods — 完整版 (V62, 2026-04-25)

## 1. Spatial Transcriptomics Analysis (Figure 1)

**Dataset**: Kuppe et al., 2022, Nature (Zenodo 6578047, Visium V1, 6 human hearts, 19,525 spots)

**Preprocessing**: 
- Visium .h5ad files processed via Python (scanpy v1.9.3) and R (Seurat v4.3.0)
- Normalization: log1p transformation of UMI counts
- Zone annotation: FZ (fibrotic zone), IZ (ischemic zone), IZ_BZ (border zone), Healthy

**Statistical Analysis**:
- Spot-level NDUFB7 positivity: spot counted as "positive" if normalized UMI > 0
- Zone comparison: Mann-Whitney U test (spot-level, n=19,525) with effect size (Cohen's d)
- **Pseudoreplication correction**: Linear mixed-effects model (lme4 v1.1-35.1) with sample ID as random effect (effective n=5 hearts)
- Spatial autocorrelation: Moran's I (spdep v1.2-8) reported for sensitivity analysis

**Key Parameters**:
- Spot diameter: 55 μm (~1-10 cells/spot)
- Positive rate threshold: normalized UMI > 0
- Random effects structure: (1|sample_id)

---

## 2. Multi-Platform Bulk Validation (Figure 2, SuppFig 3)

**Cohorts**:
| Platform | Dataset | n | Comparison | Normalization |
|----------|---------|---|------------|---------------|
| Affymetrix 1.1 ST | GSE57338 | 313 | HF vs NF | RMA (oligo v1.62.0) |
| RNA-seq RPKM | GSE116250 | 64 | DCM vs NF | RPKM (length-normalized) |
| RNA-seq counts | GSE55296 | 37 | DCM vs ICM | DESeq2 VST |
| Affymetrix U133A | GDS4772 | 17 | DCM vs Normal | RMA (affy v1.72.0) |
| Proteomics iBAQ | PXD010154 | 19 fractions | Detection | MaxQuant v1.6.0.16 |

**Statistical Framework**:
- Each platform: Wilcoxon rank-sum test (HF vs control) + Cohen's d effect size
- Cross-platform synthesis: Random-effects meta-analysis (meta v6.2-1) with I² heterogeneity
- **Length bias correction**: RPKM values flagged for 137aa length artifact; counts-based platform prioritized
- **Forest plot**: Effect sizes (Cohen's d) with 95% CI, stratified by platform

**GDS4772 Specifics**:
- Downloaded via GEOquery (v2.66.0) with local soft file caching
- GPL6244 annotation bypassed due to FTP timeout; direct GDS4772.soft parsing used
- Probe ID 8034843 mapped to NDUFB7 (verified by GDS4772_expression.txt)

---

## 3. Weighted Gene Co-expression Network Analysis (Figure 3)

**Dataset**: GSE57338 (n=313)

**Parameters** (frozen at V25):
- Soft-thresholding power: β=6 (signed network, scale-free topology R²>0.85)
- Module detection: dynamicTreeCut (deepSplit=2, minModuleSize=30)
- Module merging threshold: MEDissThres=0.25
- NDUFB7 module membership: kME=0.51 (brown module, rank 164/221)

**Validation**:
- Module stability: bootstrap resampling (n=100)
- Hub gene confirmation: top 20 kME genes in brown module

---

## 4. Expression Quantitative Trait Loci (Figure 4)

**eQTLGen Phase I** (blood, n=31,684):
- cis-window: ±1Mb from NDUFB7 TSS
- Lead SNP: rs11085898 (chr19:14641837, p=2.6×10⁻⁹)
- Effect direction: T allele increases NDUFB7 expression

**GTEx v11 Heart Left Ventricle** (n=387):
- Lead SNP: rs8103021 (chr19:14641837, p=1.89×10⁻⁶)
- Effect direction: T allele decreases NDUFB7 expression (slope=-0.080)
- **Ensembl ID**: ENSG00000099795.7 (GTEx v11 native; differs from V60-recorded ENSG00000167996)

**Cross-Tissue Comparison**:
- Directional discordance: blood ↑ vs heart ↓ for same T allele
- Colocalization: attempted with coloc v5.2.2 (PP.H4 reported where available)

---

## 5. In Silico Knockout Pathway Simulation (Figure 5)

**Method**: PROGENy (v1.20.0) on GSE57338
- Bottom 25% NDUFB7 expression (n=79) vs top 25% (n=79)
- Pathway activity: 14 PROGENy pathways (100-model bootstrap)
- Statistical test: Wilcoxon rank-sum, two-tailed

**Key Predictions**:
- OXPHOS: 9.62 vs 9.87, p=1.13×10⁻¹³
- ROS defense: 8.31 vs 8.37, p=6.06×10⁻⁴
- Fibrosis signaling: 5.93 vs 5.71, p=3.94×10⁻⁶

**Limitations**:
- In silico prediction, not experimental knockout
- HF severity may confound bottom 25% group composition
- Results labeled as "predictive" not "causal" in all figures

---

## 6. Single-Cell Reference Analysis (Supplementary)

**GSE109816** (Litviňuková et al., 2020, Nature):
- Normal human left atrium, 9,994 cells (4,987 CM + 5,007 NCM)
- UMI matrix: GSE109816_normal_heart_umi_matrix.csv.gz
- Cell annotation: GSE109816_normal_heart_cell_info.txt.gz
- **Key finding**: NCM NDUFB7 median=4 vs CM median=1 (p=4.4×10⁻⁸²)

**GSE183852** (Koenig et al., 2022, Nature Cardiovascular Research):
- Donor + DCM left ventricle, 220,752 nuclei (Seurat object)
- Condition: Donor (n=154,083) + DCM (n=66,669)
- Cell types: 15 Names annotations (Cardiomyocytes, Fibroblasts, Endothelium, etc.)
- **Key finding**: NDUFB7 extremely sparse (87.9% zero-value, median=0, mean=0.14)
- **Anatomical caveat**: Atrium (GSE109816) vs ventricle (GSE183852) may confound comparison

**Cross-Dataset Comparison**:
- Normalization: log1p(UMI+1) for visualization
- Statistical test: Wilcoxon rank-sum (conservative, non-parametric)
- Cell type mapping: published Names annotations (GSE183852) vs N_LA_CM/N_LA_NCM (GSE109816)

---

## 7. Deconvolution Analysis (V61, SuppFig 6-8)

**Method**: Published cell2location deconvolution (Kuppe et al., 2022, Nature)
- Reference: Kuppe 2022 scRNA-seq (CM 28.9%, FB 7.9%, EC 15.5%, MP 8.1%)
- Algorithm: cell2location v0.1.3 (PyTorch-based)
- Spot-level cell type proportions: 11 cell types

**Key Results**:
- NDUFB7-high spots: CM enrichment +5.6% (mean, 4/5 samples)
- NDUFB7-high spots: FB depletion -3.7%
- Spearman correlation (CM vs NDUFB7): mean rho=0.226, p<0.001

---

## 8. Software and Code

**R environment**: v4.3.1
**Key packages**: 
- Seurat v4.3.0, SeuratObject v4.1.3
- ggplot2 v3.4.2, data.table v1.14.8
- lme4 v1.1-35.1, spdep v1.2-8
- GEOquery v2.66.0, Biobase v2.58.0
- PROGENy v1.20.0, coloc v5.2.2

**Python environment**: v3.10.9
**Key packages**: scanpy v1.9.3, pandas v2.0.3, numpy v1.24.3

**Code availability**: 
- GitHub: https://github.com/sots-melb/NDUFB7-HF-2026
- Zenodo DOI: [pending]
- License: MIT

---

## 9. Statistical Transparency

**p-value reporting**: All p-values reported with:
- Test type (Wilcoxon/Mann-Whitney/LMM)
- Sample size (n)
- Effect size (Cohen's d or beta)
- Raw p-value
- **FDR-corrected q-value** (BH method, where applicable)

**Missing data**: 
- GSE183852 cell type-specific statistics pending full matrix extraction (Robj memory constraints)
- SMR/HEIDI analysis pending HERMES GWAS outcome data availability

---

*Methods V62 — 2026-04-25*
*Frozen for Revision preparation*
