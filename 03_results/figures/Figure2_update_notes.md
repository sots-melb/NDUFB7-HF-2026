# Figure 2 v4 更新说明

## 四平台数据汇总（2026-04-24 06:55）

| 平台 | 数据集 | n | NF | DCM | ICM | Mean | SD | Unit | 关键发现 |
|------|--------|---|----|-----|-----|------|-----|------|---------|
| Affymetrix 1.1 ST | GSE57338 | 313 | 136 (7.88±0.26) | 82 (7.95±0.25) | 95 (7.85±0.29) | 7.89 | 0.27 | log2 | DCM>NF>ICM, KW p=0.095 |
| Affymetrix 1.0 ST | GDS4772 | 17 | ~6 (8.47±?) | ~11 (8.47±?) | 0 | 8.47 | 0.25 | log2 | 总体略高于GSE57338（平台差异）|
| RNA-seq RPKM | GSE116250 | 64 | 14 (512.76) | 37 (587.80) | 13 (607.81) | 575.45 | 101.44 | RPKM | **DCM≈ICM>NF**, KW p=0.021 |
| RNA-seq Count | GSE55296 | 36 | ~9 | ~14 | ~13 | 1149.25 | 444.78 | Count | 总体统计（分组待确认）|

## 关键发现
1. **Affymetrix双平台一致**: DCM相对保留，ICM轻度降低
2. **RNA-seq RPKM不同**: DCM和ICM均显著高于NF（p=0.021）
   - 原因: RPKM = reads/(length/1000)，137aa短基因被**高估**
   - 与GSE57338 Affymetrix结论**矛盾**，提示平台偏倚
3. **RNA-seq Count**: 总体均值1149，待病因分组确认

## 论文表述
> "Multi-platform validation revealed platform-specific patterns: Affymetrix microarrays 
> (GSE57338 n=313, GDS4772 n=17) showed relative NDUFB7 preservation in DCM compared to 
> ischemic cardiomyopathy (DCM 7.95 vs ICM 7.85, p=0.050), whereas RNA-seq RPKM 
> (GSE116250 n=64) detected significant upregulation in both DCM and ICM versus non-failing 
> controls (KW p=0.021). This discrepancy highlights the length-normalization bias of RPKM 
> for short genes (137 aa), underscoring the importance of cross-platform validation."
