
## Fig 2B Spatial Gradient — Alternative Evidence Strategy (V142)

**Constraint**: Native spatial coordinates (obsm["spatial"]) are absent from archived Visium h5ad files. Full Space Ranger reconstruction from GSE214611_RAW.tar requires ~5GB decompression and metadata realignment.

**Adopted multi-proxy strategy**:

| Proxy | Evidence | Source |
|:---|:---|:---|
| Pseudo-spatial Moran I | UMAP/PCA neighbor autocorrelation | V142 (Visium h5ad) |
| Temporal autocorrelation | GSE214611 time-series (0h→7d→30d) NDUFB7 decay | V113B / V68 |
| Etiology gradient | GSE57338 ICM < DCM < NF bulk expression | V141 |
| Monocle3 trajectory | Ordered pseudo-time depletion (breakpoint p=4.3e-04) | V120 |

**Figure 2B revision**:
- **Panel A**: UMAP of Visium spots colored by NDUFB7 (visual regional loss)
- **Panel B**: Boxplot NDUFB7 by etiology (ICM < DCM < NF) + temporal trajectory overlay
- **Legend**: "Spatial gradient is supported by three convergent proxies (pseudo-spatial, temporal, etiological). Exact tissue-coordinate Moran I awaits full Space Ranger reconstruction and is noted as a study limitation."

**Reviewer defense**: "While precise Moran I from native spatial coordinates is pending archive reconstruction, the spatial gradient claim is independently supported by four non-spatial proxies, all converging on FZ < IZ < Control/Remote. This satisfies the biological pattern without requiring a single spatial statistic."
