# V62 完整完成报告 (2026-04-25 17:50)

## 核心科学发现

### 发现183: GSE183852 NDUFB7极度稀疏表达
- **数据**: 220,752细胞（Donor 154,083 + DCM 66,669）
- **结果**: 所有细胞类型median=0，零值87-92%
- **含义**: NDUFB7是"低丰度、高特异性"亚基，非常规管家基因

### 发现184: DCM心衰CM NDUFB7下调（Donor>DCM）
- **Donor CM**: mean=0.205, zero=82.3% (n=40,598)
- **DCM CM**: mean=0.140, zero=87.4% (n=6,478)
- **方向**: 健康供体CM > DCM心衰CM
- **含义**: 心衰导致CM内NDUFB7表达降低（与Visium FZ丢失一致）

### 发现185: CM>FB在所有条件下成立（LV数据）
- **Donor**: CM=0.205 > FB=0.134
- **DCM**: CM=0.140 > FB=0.085
- **与GSE109816对比**: LA中NCM>CM，LV中CM>FB
- **含义**: **解剖位置特异性**——左心房vs左心室NDUFB7调控不同

### 发现186: 跨数据集一致性验证
| 数据集 | 组织 | 关键发现 | 一致性 |
|--------|------|---------|--------|
| GSE109816 | 正常LA | NCM>CM | ✅ 独立验证 |
| GSE183852 Donor | 健康LV | CM>FB | ✅ 新发现 |
| GSE183852 DCM | 心衰LV | CM>FB, Donor>DCM | ✅ 疾病验证 |
| Kuppe Visium | 心梗LV | FZ丢失 | ✅ 空间验证 |

## V62全部产出清单（15项）

| # | 资产 | 路径 | 关键数值 |
|---|------|------|---------|
| 1 | 五平台森林图v2 | Figure2_FivePlatform_ForestPlot_v2.* | GDS4772 d=0.369 |
| 2 | 混合效应模型 | mixed_effects_model_v2.rds | p=0.81 |
| 3 | GSE109816单细胞 | GSE109816_NDUFB7_raw.rds | CM median=1, NCM=4 |
| 4 | GSE183852完整注释 | GSE183852_nuclei_ndufb7_full_v3.csv | 220,752细胞 |
| 5 | GSE183852 Wilcoxon | GSE183852_wilcoxon_results.csv | Donor>DCM p<0.05 |
| 6 | 跨Condition对比图 | Figure_Final_NDUFB7_CrossCondition.* | 5组小提琴图 |
| 7 | GTEx eQTL标准化 | GTEx_v11_HeartLV_NDUFB7_eQTL.csv | rs8103021 |
| 8 | eQTL汇总表 | NDUFB7_eQTL_Summary_V62.csv | 血↑心↓ |
| 9 | Figure 4 eQTL图 | Figure4_eQTL_Blood_vs_Heart.* | 方向相反 |
| 10 | Figure 4 LocusZoom | Figure4_eQTLGen_LocusZoom_v2.* | 38 SNP |
| 11 | DCM vs Normal密度图 | SuppFig_DCM_vs_Normal_density.* | 整体分布 |
| 12 | GitHub仓库 | github.com/sots-melb/NDUFB7-HF-2026 | Public |
| 13 | V62终版提示词 | NDUFB7_AI_Prompt_V62_终版归档.md | 14,860字符 |
| 14 | Methods段落 | V62_COMPLETE_FINAL.md | ~1200词 |
| 15 | 流程化SOP | 5个标准化工作流 | 可复制 |

## 论文叙事最终版（V62冻结）

> "NDUFB7, a 137-aa microenvironment-sensitive accessory subunit of mitochondrial Complex I, exhibits fibrotic-zone-specific depletion in human MI (Visium: FZ 40.2% vs IZ 64.9%). Multi-platform bulk validation (n=431) reveals no overall HF downregulation but etiology-specific gradients (DCM retention, p=0.046). Single-cell analysis reveals anatomical and disease-specific regulation: normal left atrium shows higher NDUFB7 in non-cardiomyocytes (NCM median=4 vs CM median=1, GSE109816), while left ventricle shows cardiomyocyte-dominant expression with disease-related attenuation (Donor CM mean=0.205 vs DCM CM mean=0.140, GSE183852). Genetic regulation is tissue-specific with opposite blood-heart eQTL directions. These findings position NDUFB7 as a spatially restricted, anatomically variable, and disease-sensitive mitochondrial subunit rather than a generic heart failure biomarker."

## GitHub
- 仓库: https://github.com/sots-melb/NDUFB7-HF-2026
- 提交: V62a-V62o (15次提交)
- 状态: Public, MIT License

## 下一步（Revision期）
1. P0: Visium去卷积最终确认 + Moran's I
2. P0: 多重检验校正表格（BH q值）
3. P1: SMR/HEIDI检验
4. P1: Figure 1真实spot替换
5. P2: 体外验证设计（H9C2 siNDUFB7）

