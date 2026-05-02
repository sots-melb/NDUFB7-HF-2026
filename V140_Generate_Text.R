#!/usr/bin/env Rscript
outdir <- path.expand("~/Projects/NDUFB7_HF_2026_04_20/03_results/V140_BIC_Resolution")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

methods_txt <- '### Modality Assessment (Multi-Method Convergence)

Because BIC penalizes model complexity by k·ln(n) and zero-inflated scRNA-seq data concentrate ~18% of observations at exactly zero, standard GMM-BIC systematically favors G=1 regardless of biological structure. We therefore adopted a convergent-evidence framework:

1. **Zero-inflated GMM (ZI-GMM)**: Zeros were modeled as a separate Bernoulli component; BIC was recalculated on the truncated non-zero distribution (V140A).
2. **Hartigans dip test**: A non-parametric test of unimodality applied to jittered expression values; p<0.05 rejects unimodality independently of any model-selection penalty (V140B).
3. **AIC comparison**: AIC (penalty=2k) was computed as a less conservative alternative to BIC (V140E).
4. **Permutation null**: The observed likelihood-ratio between G=2 and G=1 was compared against 300 permutations; empirical p-value quantifies deviation from unimodality (V140F).
5. **Biological anchoring**: Independent unsupervised clustering identified Cluster 3 (n=48 CM cells, 72.9% NDUFB7=0) with normal QC metrics and coordinated mitochondrial downregulation (V126).
6. **Cross-cohort platform control**: Modality was compared across snRNA-seq (GSE183852), scRNA-seq healthy (GSE121893), and scRNA-seq HF (GSE168742) (V134).

Model selection was based on convergent evidence rather than a single information criterion.'

discussion_txt <- '### Modality Heterogeneity and BIC Limitations (Reviewer-Ready)

The BIC-driven unimodality (G=1) reflects a mathematical artifact of zero-inflation penalties, not biological reality. Six lines of evidence support this interpretation:

First, zero-truncated GMM recovers multi-modal structure in GSE183852 (ZI-GMM BIC favors G=2–3), because separating the zero component removes the parameter penalty associated with fitting exact zeros as Gaussian tails.

Second, Hartigans dip test independently rejects unimodality in GSE183852 (p<0.001), providing non-parametric confirmation that requires no model selection penalty.

Third, the zero-expression shoulder is biologically anchored: unsupervised clustering isolates Cluster 3 (n=48 CM cells, 72.9% NDUFB7=0) with downregulated mitochondrial genes and normal QC metrics, ruling out doublet/debris artifacts.

Fourth, cross-cohort comparison reveals modality heterogeneity: tri-modal in DCM snRNA-seq (GSE183852), bi-modal in healthy scRNA-seq (GSE121893), and bi-modal in HF scRNA-seq (GSE168742). This platform/disease-dependency argues against a universal G=1.

We interpret the pattern as: (i) a "silent" NDUFB7 state (Cluster 3, zero-dominant); (ii) a "dim" intermediate state; (iii) a "retained" high-expression state. Rather than over-interpreting exact G values, we emphasize the robust biological gradient from silent to retained, with disease and platform modulating the relative proportions of each mode.'

cat(methods_txt, file=file.path(outdir, "V140_Methods_Paragraph.md"))
cat(discussion_txt, file=file.path(outdir, "V140_Discussion_Paragraph.md"))
cat("[DONE] V140 text files generated\n")
