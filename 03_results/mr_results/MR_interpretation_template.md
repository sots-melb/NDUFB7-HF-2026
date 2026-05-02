
# NDUFB7 MR分析结果解读模板

## 暴露数据
- 数据源: eQTLGen Phase I (31,684 samples, 16 cohorts)
- 基因: NDUFB7 (ENSG00000099795, chr19:14679882)
- cis-eQTL SNP数: 38 (F>10)
- 最强SNP: rs11085898 (p=2.64e-9, Z=5.95)

## 结局数据
- 待填入: IEU GWAS ID
- 表型: Heart failure
- 样本量: 待填入

## MR结果（待填入）
| Method | SNPs | Beta | SE | P-value | 解读 |
|--------|------|------|-----|---------|------|
| Wald Ratio | ? | ? | ? | ? | 单SNP因果估计 |
| IVW | ? | ? | ? | ? | 主分析（多SNP加权） |
| MR-Egger | ? | ? | ? | ? | 水平多效性校正 |
| Weighted Median | ? | ? | ? | ? | 50%有效IV假设 |

## 敏感性（待填入）
- 异质性Q_pval: ?
- MR-Egger intercept p: ?
- 结论: ?

## 与论文整合
如果MR显著（p<0.05）:
> "Mendelian randomization using 38 cis-eQTL instruments (F>10) revealed a 
> [protective/risk-increasing] causal association between genetically predicted 
> NDUFB7 expression and heart failure (IVW beta=?, p=?), supporting [机制方向]."

如果MR不显著:
> "MR analysis did not detect a significant causal effect (IVW p=?), suggesting 
> that the observed spatial redistribution may be a consequence rather than cause 
> of heart failure, or that the effect is context-dependent (cell-type specific)."
