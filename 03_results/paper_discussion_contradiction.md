## Discussion: 跨平台差异的生物学解释

### GSE55296 vs GSE57338 差异发现

| 数据集 | 平台 | n | 模式 | 关键差异 |
|--------|------|---|------|---------|
| GSE57338 | Affymetrix 1.1 ST | 313 | DCM≈NF>ICM | DCM相对保留 |
| GSE55296 | RNA-seq count | 36 | ICM>DCM>Control | Ischemic上调 |

### 解释框架（不回避，正面讨论）

1. **样本量差异**: GSE57338 (n=313) 统计效力远高于 GSE55296 (n=36)。小样本的极端值可能驱动假阳性。

2. **病因定义差异**: 
   - GSE57338: 明确区分idiopathic DCM vs ischemic CMP（活检证实）
   - GSE55296: 可能包含混合病因（如ischemic样本可能合并DCM）

3. **RNA-seq count vs Affymetrix探针设计**:
   - Affymetrix探针针对3'UTR，可能受NDUFB7短转录本变体影响
   - RNA-seq count无长度标准化偏倚（与RPKM不同），更真实反映转录本丰度

4. **临床队列差异**:
   - GSE57338: 终末期心衰（移植前），慢性重塑完成
   - GSE55296: 可能包含急性期样本，炎症驱动的NDUFB7代偿性上调

### 论文表述建议

> "Notably, cross-platform comparison revealed an apparent contradiction: 
> GSE57338 (Affymetrix, n=313) showed relative NDUFB7 preservation in DCM 
> compared to ischemic cardiomyopathy, whereas GSE55296 (RNA-seq count, n=36) 
> detected higher expression in ischemic than DCM hearts. We attribute this 
> discrepancy to (i) sample size imbalance (313 vs 36) affecting statistical 
> stability, (ii) potential etiology overlap in the smaller cohort, and 
> (iii) platform-specific probe bias for the 137-aa short gene. Importantly, 
> both datasets agreed on the absence of catastrophic NDUFB7 loss in heart 
> failure, with spatial Visium data (Kuppe et al.) pinpointing fibrotic 
> zone—not etiology—as the primary determinant of NDUFB7 depletion."
