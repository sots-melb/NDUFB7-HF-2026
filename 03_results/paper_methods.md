# Methods

## Spatial Transcriptomics Analysis

**Dataset**: Kuppe et al. Nature 2022 human MI Visium dataset (Zenodo 6578047, cellxgene.cziscience.com). Six samples: 3 post-MI (P3, P7, P14, P20) and 3 healthy controls (P5, P13).

**Processing**: .h5ad files were loaded via Scanpy (v1.9). NDUFB7 expression was extracted using Ensembl ID ENSG00000099795. Spot-level detection was defined as normalized expression > 0. Cross-zone comparisons used Mann-Whitney U test (scipy.stats.mannwhitneyu).

## Bulk Transcriptomic Validation

**Datasets**: 
- GSE57338 (Affymetrix Human Gene 1.1 ST, n=313 left ventricle samples: 136 non-failing, 82 DCM, 95 ICM)
- GDS4772 (Affymetrix Human Gene 1.0 ST, n=17: 5 normal, 12 DCM)
- GSE116250 (RNA-seq RPKM, n=64: 14 NF, 37 DCM, 13 ICM)
- GSE55296 (RNA-seq count, n=36: 10 control, 13 DCM, 13 ischemic)

**Analysis**: NDUFB7 was extracted using Ensembl ID ENSG00000099795 or gene symbol matching. Group comparisons used Kruskal-Wallis test (kruskal.test) and pairwise Wilcoxon rank-sum test (wilcox.test) in R v4.3. Effect sizes reported as Cohen's d or rank-biserial correlation.

## Complex I Network Analysis

From GSE57338, all 11 B-subfamily members (NDUFB1-11) plus 4 reference subunits (NDUFS4, NDUFA9, NDUFV1, NDUFC2) were extracted via GPL11532 platform annotation. Pearson correlation was computed across 313 samples. Hierarchical clustering and heatmap visualization used ggplot2.

## Genetic Regulation Analysis

**eQTLGen Phase I**: cis-eQTL summary statistics (2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz) were downloaded from eqtlgen.org. NDUFB7 (ENSG00000099795) cis-eQTLs were extracted by grep. F-statistics computed as (beta/se)².

**GTEx v11**: Heart Left Ventricle eQTLs were downloaded from GTEx Portal (GTEx_Analysis_v11_eQTL.tar). eGenes and significant variant-gene pairs extracted using zcat/grep and pyarrow (parquet format).

## Protein Validation

PXD010154 human cardiac proteome atlas (PRIDE) was searched for NDUFB7 peptides. Coverage calculated as unique peptides / total amino acids (137 aa).

## Statistical Standards

Significance threshold: p < 0.05. Multiple testing correction: FDR (Benjamini-Hochberg) where applicable. All analyses performed in R v4.3.1 or Python 3.12.

## Data Availability

All data are publicly available via GEO, PRIDE, GTEx Portal, and eQTLGen consortium. Analysis scripts are available upon request.
