## Detailed Methods Parameters

### 2.1 Data Preprocessing
- **GSE57338 (GPL570):** RMA normalization using affy package v1.76.0; batch correction by ComBat v3.42.0 (sva package); log2 transformation.
- **GSE59867 (GPL11532):** RMA normalization using oligo package v1.62.0; fRMA batch correction; probe-to-gene mapping using Brainarray custom CDF v25.0.0.
- **GSE168742:** Cell Ranger v7.1.0 alignment to mm10; Seurat v4.3.0 QC (nFeature_RNA 200-8000, percent.mt <15%); SCTransform normalization; Harmony v0.1.1 batch integration.
- **GSE157282 / GSE243655:** featureCounts v2.0.1 against GENCODE v34; TPM normalization using edgeR v3.40.0 calcNormFactors + limma v3.54.0 voom.

### 2.2 Statistical Thresholds
- **Genome-wide significance:** p < 5×10⁻⁸ (GWAS); p < 0.05/n_tests (Bonferroni); FDR < 0.05 (Benjamini-Hochberg).
- **Targeted hypothesis tests:** Nominal p < 0.05 (pre-specified NDUFB7-death relationships).
- **Effect size interpretation:** Cohen's d: small 0.2, medium 0.5, large 0.8; Spearman rho: weak <0.3, moderate 0.3-0.5, strong >0.5.

### 2.3 SMR Parameters
- **Software:** SMR v1.3.1 Linux x86_64.
- **eQTL source:** GTEx v8 Heart-LV (n=372); eQTLGen blood (n=31,456).
- **GWAS source:** HERMES GWAS (n=977,323; 47 cohorts).
- **Thresholds:** SMR p < 0.05; HEIDI p > 0.05 (excludes LD confounding); minimum 3 SNPs in HEIDI test.
- **Instrument selection:** Top eQTL SNP per gene (lowest p-value); F-statistic > 10 for instrument strength.

### 2.4 Partial Correlation
- **ROS control genes:** NOX1, NOX4, SOD1, SOD2, CAT, PRDX1, TXNRD1, GCLC, GCLM (n=9).
- **Method:** ppcor package v1.1; Spearman partial correlation controlling for ROS gene expression matrix.
- **Attenuation metric:** delta = |rho_raw| - |rho_partial|; attenuation > 0.05 considered ROS-mediated.

### 2.5 Spatial Analysis
- **Visium preprocessing:** Space Ranger v2.0.0; SpotQC (nUMI > 500, nGene > 250).
- **Normalization:** SCTransform (v2 regularization); integration by Harmony.
- **Regional annotation:** Manual annotation based on H&E morphology + marker genes (ACTA2 for FZ, TNNI3 for viable CM).
- **Spatial statistics:** Moran's I using spdep v1.2-8; Kruskal-Wallis for regional comparisons.

### 2.6 Software Environment
- **OS:** Ubuntu 22.04 LTS (x86_64).
- **R:** v4.3.1 (R Foundation).
- **Python:** v3.10.12.
- **Key packages:** data.table_1.14.8, dplyr_1.1.3, ggplot2_3.4.3, Seurat_4.3.0, smr_1.3.1.
