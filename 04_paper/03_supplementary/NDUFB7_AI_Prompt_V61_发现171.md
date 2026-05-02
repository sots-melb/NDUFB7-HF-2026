# V61增量: 发现171 — 混合效应模型纠正伪重复 (2026-04-25 08:43)

## 发现171详情

| 属性 | 内容 |
|------|------|
| 名称 | 混合效应模型纠正Visium伪重复 |
| 原始结论 | FZ 40.2% vs IZ 64.9%, p<0.0001 (Mann-Whitney, n=16004 spots) |
| 修正结论 | zoneFZ p=0.805, zoneIZ p=0.950 (LMM, sample随机效应, df≈5) |
| 效应量 | FZ vs Healthy Cohen's d = -0.243 |
| 科学等级 | 从Tier 1确证 → 降级为Tier 2初步证据 |

## 各组统计 (修正版)

| 区域 | 样本 | mean | median | n_spots |
|------|------|------|--------|---------|
| Healthy | GT_P13/P15/P9 | 1.098 | 1.289 | 9823 |
| IZ | P3 | 1.043 | 1.260 | 3771 |
| FZ | P20 | 0.858 | 0.000 | 2410 |

## 核心洞察

1. FZ的median=0是真实观察 (2410 spots中大量NDUFB7=0)
2. 但样本间变异(随机效应方差=0.208)解释了组间差异
3. 当前5样本不足以支持"纤维化区特异性"的统计推断
4. **Reviewer #1 Major 1质疑完全正确**

## 叙事调整 (必须执行)

### 原表述 (禁止)
"NDUFB7 exhibits fibrotic-zone-specific depletion in human MI (Visium: FZ 40.2% vs IZ 64.9%, p<0.0001)"

### 修正表述 (必须)
"Spot-level analysis suggested lower NDUFB7 detection in the fibrotic zone (median=0.00 vs infarct zone median=1.26), but this pattern was not statistically significant after accounting for inter-sample clustering via linear mixed-effects modeling (β=-0.148, p=0.81, n=5 hearts). This positions NDUFB7 spatial redistribution as a hypothesis requiring validation in larger cohorts."

## 对Reviewer回应的影响

Reviewer #1 Major 1 (原): "Visium伪重复问题"
回应: ✅ 完全接受质疑，LMM纠正后差异消失，诚实报告

Reviewer #1 Major 2 (原): "四平台矛盾"
回应: 现在更一致了——Bulk无差异 + 空间不显著 = "NDUFB7不是强HF标志物"

## 项目影响评估

| 维度 | 影响 |
|------|------|
| 核心叙事 | 从"空间丢失"转向"病因特异性+方法学警示" |
| 证据等级 | 空间发现从Tier 1 → Tier 2 |
| 论文诚实性 | ⬆️ 大幅提升 |
| 接收概率 | ⬆️ 审稿人更欣赏自我纠正 |
| 国自然标书 | 需调整假说：从"空间标志物"→"病因敏感亚基" |

## 下一步任务调整

🔴 P0 (本周):
1. GDS4772完成 → 五平台森林图 (Reviewer #1 Major 2)
2. 论文叙事全面调整 (空间发现降级)
3. TwoSampleMR本地VCF或Coloc (替代SMR)

🟡 P1 (下周):
4. GSE109816剪接探索
5. 代码GitHub上传

🟢 P2:
6. 更大样本空间验证 (未来队列合作)

