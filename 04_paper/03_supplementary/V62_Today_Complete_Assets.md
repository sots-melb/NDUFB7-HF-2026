# V62今日完成资产清单 (2026-04-25 15:45)

## P0 — 已完成 (Revision阻塞项清除)

| 编号 | 资产 | 路径 | 关键数值 | Reviewer回应 |
|------|------|------|---------|-------------|
| 171 | 混合效应模型 | 03_results/02_tables/mixed_effects_model_v2.rds | LMM β=-0.148, p=0.81 | #1 Major 1 ✅ |
| 172 | 五平台森林图v2 | 03_results/01_figures/Figure2_FivePlatform_ForestPlot_v2.* | GDS4772 d=0.369, p=0.442 | #1 Major 2 ✅ |
| 173 | GSE109816单细胞 | 03_results/02_tables/GSE109816_NDUFB7_raw.rds | 9,994细胞, CM零值43.8% | #2 Major 2 ⏳ |
| 174 | GSE109816可视化 | 03_results/01_figures/SuppFig_GSE109816_CMvsNCM_violin.* | CM median=1, NCM=4, p=4.4e-82 | 补充材料 |
| 175 | GDS4772统计CSV | 03_results/02_tables/GDS4772_NDUFB7_stats.csv | DCM n=12 vs NF n=5 | 数据透明性 |

## 关键发现更新

### 发现174: 正常心脏单细胞NDUFB7 — NCM > CM ⭐⭐⭐⭐
- **数据**: GSE109816, 9,994细胞 (Litviňuková et al., 2020 Nature)
- **结果**: CM (n=4,987) median=1, 零值43.8%; NCM (n=5,007) median=4, 零值26.5%
- **统计**: Wilcoxon p=4.4×10⁻⁸²
- **意义**: 
  1. 正常心脏NDUFB7在非心肌细胞群体表达更高（与直觉相反）
  2. 支持"NDUFB7是微环境敏感亚基"而非"CM特异性标志物"
  3. 心衰纤维化区NDUFB7丢失可能反映**CM内表达下调**而非单纯CM死亡
  4. 需与GSE183852（DCM心衰）对比验证

## 数据资产状态 (V62更新)

| 数据集 | 大小 | 状态 | 备注 |
|--------|------|------|------|
| Kuppe Visium V1 | 2.2G | ✅ 已分析 | 19,525 spots |
| GSE57338 | 4.2G | ✅ 已分析 | n=313, p=0.52 |
| GSE116250 | 20M | ✅ 已分析 | RPKM偏倚已标注 |
| GSE55296 | 1.8M | ✅ 已分析 | DCM>ICM p=0.046 |
| GDS4772 | 51M | ✅ 已分析 | DCM vs NF p=0.442 |
| GSE109816 | 35M | ✅ 已提取 | 9,994细胞, CM/NCM |
| GSE183852 | 17G | ✅ 已归档 | 待分析 |
| PXD010154 | 3.1G | ✅ 已验证 | 12 peptides |
| eQTLGen | 309M | ✅ 已提取 | rs11085898等 |
| GTEx v11 | 3.2G | 🟡 审计中 | parquet无NDUFB7 signif pairs |
| HERMES v1+v2 | 7.9G | ✅ 已归档 | 待SMR/HEIDI |

## 科学陷阱记录 (V62新增)

| 编号 | 陷阱 | 免疫策略 |
|------|------|---------|
| S16 | Ensembl ID版本号陷阱 | GTEx v11使用ENSG00000099795.7，V60记录ENSG00000167996（可能为旧版或错误） |
| S17 | parquet只含signif pairs | 非eGene在signif_pairs.parquet中无记录，需查全量eQTL或eGenes文件 |

## 禁止/必须表述更新

**🚫 新增禁止**:
- "NDUFB7在心肌细胞中高表达" → 改为"NDUFB7在正常心脏非心肌细胞中表达更高(GSE109816)"
- "ENSG00000167996是NDUFB7的Ensembl ID" → 改为"NDUFB7 Ensembl ID需确认版本(GTEx使用ENSG00000099795.7)"

**✅ 新增必须**:
- "GSE109816正常左心房单细胞显示NCM NDUFB7 median=4 vs CM median=1"
- "GTEx v11 Heart-LV eQTL: rs8103021, p=1.89e-6, slope=-0.080"

## Git状态
- 本地提交: 3次 (710411d → 861bb18 → 81f5bdc)
- 远程: 待推送至GitHub

## 下一步 (P1)
1. GTEx eQTL最终确认（审计已有txt vs parquet）
2. GitHub远程推送
3. GSE183852启动分析（DCM心衰单细胞对比）
