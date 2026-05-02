### Modality Assessment (Multi-Method Convergence)

Because BIC penalizes model complexity by k·ln(n) and zero-inflated scRNA-seq data concentrate ~18% of observations at exactly zero, standard GMM-BIC systematically favors G=1 regardless of biological structure. We therefore adopted a convergent-evidence framework:

1. **Zero-inflated GMM (ZI-GMM)**: Zeros were modeled as a separate Bernoulli component; BIC was recalculated on the truncated non-zero distribution (V140A).
2. **Hartigans dip test**: A non-parametric test of unimodality applied to jittered expression values; p<0.05 rejects unimodality independently of any model-selection penalty (V140B).
3. **AIC comparison**: AIC (penalty=2k) was computed as a less conservative alternative to BIC (V140E).
4. **Permutation null**: The observed likelihood-ratio between G=2 and G=1 was compared against 300 permutations; empirical p-value quantifies deviation from unimodality (V140F).
5. **Biological anchoring**: Independent unsupervised clustering identified Cluster 3 (n=48 CM cells, 72.9% NDUFB7=0) with normal QC metrics and coordinated mitochondrial downregulation (V126).
6. **Cross-cohort platform control**: Modality was compared across snRNA-seq (GSE183852), scRNA-seq healthy (GSE121893), and scRNA-seq HF (GSE168742) (V134).

Model selection was based on convergent evidence rather than a single information criterion.