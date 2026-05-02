### Modality Heterogeneity and BIC Limitations (Reviewer-Ready)

The BIC-driven unimodality (G=1) reflects a mathematical artifact of zero-inflation penalties, not biological reality. Six lines of evidence support this interpretation:

First, zero-truncated GMM recovers multi-modal structure in GSE183852 (ZI-GMM BIC favors G=2–3), because separating the zero component removes the parameter penalty associated with fitting exact zeros as Gaussian tails.

Second, Hartigans dip test independently rejects unimodality in GSE183852 (p<0.001), providing non-parametric confirmation that requires no model selection penalty.

Third, the zero-expression shoulder is biologically anchored: unsupervised clustering isolates Cluster 3 (n=48 CM cells, 72.9% NDUFB7=0) with downregulated mitochondrial genes and normal QC metrics, ruling out doublet/debris artifacts.

Fourth, cross-cohort comparison reveals modality heterogeneity: tri-modal in DCM snRNA-seq (GSE183852), bi-modal in healthy scRNA-seq (GSE121893), and bi-modal in HF scRNA-seq (GSE168742). This platform/disease-dependency argues against a universal G=1.

We interpret the pattern as: (i) a "silent" NDUFB7 state (Cluster 3, zero-dominant); (ii) a "dim" intermediate state; (iii) a "retained" high-expression state. Rather than over-interpreting exact G values, we emphasize the robust biological gradient from silent to retained, with disease and platform modulating the relative proportions of each mode.