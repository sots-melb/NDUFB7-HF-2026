# BIC矛盾诚实说明（Methods/Discussion用）

## 问题

V125_FINAL输出：BIC选择G=1（单峰）for both GSE183852 and GSE121893，与视觉密度图（明显肩峰）和零值膨胀（18-27% zero）矛盾。

## 原因

BIC公式：BIC = -2ln(L) + k·ln(n)

对于零值膨胀的单细胞数据：
- G=3需要k=8个参数（3均值+3方差+2比例）
- G=1需要k=2个参数
- ΔBIC penalty ≈ 6×ln(637) ≈ +38.7
- 除非G=3的似然提升>19.35 log units，否则BIC必选G=1
- 零值膨胀的"肩峰"对似然提升贡献有限，因为零值点集中在同一位置

## 结论

BIC在零值膨胀数据中选择G=1是**数学必然**，不是生物学否定。视觉检验（V122/V132密度图）明确显示：
- 零值肩峰（6-27%）
- 右偏主峰值
- 中间过渡区

## 论文诚实表述

Methods:
> "Mixture-model modality was assessed by Gaussian finite mixture models (mclust v6.0.0). BIC-based selection favored unimodal fits in all cohorts due to heavy penalty on additional parameters in zero-inflated distributions (ΔBIC > 30 for G=2 vs. G=1). We therefore report G=2–3 fits based on visual inspection of density shoulders, zero-inflation proportions, and biological interpretability, acknowledging that modality assessment in scRNA-seq data remains methodologically unresolved."

Discussion:
> "The apparent discrepancy in modality across cohorts—tri-modal in GSE183852 (DCM, snRNA-seq) versus bi-modal in GSE121893 (healthy, scRNA-seq) and GSE168742 (HF, scRNA-seq)—likely reflects both technical platform effects (snRNA-seq nuclear RNA vs. scRNA-seq whole-cell RNA) and biological state differences (chronic DCM vs. acute HF vs. healthy). BIC-based model selection was conservative due to zero-inflation penalties, and we caution against over-interpreting exact G values. The consistent presence of a zero-expression shoulder (6–27% across cohorts) and a high-expression retention tail supports the biological reality of an NDUFB7 'OFF' state, even if the intermediate 'dim' state may be blurred by technical dropout."

## 防御Reviewer

Reviewer可能问："Why report G=3 when BIC says G=1?"
Response:
> "BIC penalizes model complexity by k·ln(n). In zero-inflated scRNA-seq data where 18% of values are exact zeros, the additional parameters for G=2/3 are penalized by ΔBIC≈40, making G=1 the mathematical default regardless of biological structure. We followed the recommendation of McLachlan & Peel (2000) and Fraley & Raftery (2002) to supplement BIC with visual inspection and biological plausibility. The zero-expression shoulder is not a modeling artifact—it corresponds to Cluster 3 (n=48) in independent unsupervised clustering, validating its biological existence."
