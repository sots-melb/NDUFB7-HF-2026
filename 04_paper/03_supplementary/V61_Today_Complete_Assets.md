# V61今日完成资产清单 (2026-04-25)

## 核心分析完成

| 编号 | 资产 | 路径 | 用途 | Reviewer回应 |
|------|------|------|------|-------------|
| 171 | 混合效应模型 | 03_results/02_tables/mixed_effects_model_v2.rds | LMM纠正伪重复 | #1 Major 1 ✅ |
| 172 | 健康心脏异质性 | V61提示词归档 | 10-64%零值变异 | 方法学警示 |
| 173 | 样本设计缺陷 | V61提示词归档 | n=1 FZ无法推断 | 未来实验设计 |
| - | 四/五平台森林图 | 03_results/01_figures/Figure2_FivePlatform_ForestPlot.* | 跨平台一致性 | #1 Major 2 ⏳ |
| - | TwoSampleMR IV | 03_results/02_tables/NDUFB7_IV_candidates.rds | rs7508201 | #2 Major 1 ⏳ |
| - | HERMES GWAS归档 | 01_data/03_mr_data/gwas_summary/ | v1+v2+IEU VCF | SMR/HEIDI/MR |
| - | eQTLGen Top 10 | 03_results/mr_results/ | rs11085898等 | 遗传学工具变量 |
| - | GTEx Heart-LV 6 SNP | 03_results/mr_results/ | rs8103021等 | 组织特异性eQTL |

## 数据资产状态

| 数据集 | 路径 | 大小 | 状态 |
|--------|------|------|------|
| Kuppe Visium V1 | 01_data/04_spatial_geo/Kuppe_Nature_2022/ | 2.2G | ✅ 已分析 |
| GSE57338 | 01_data/01_raw_geo/GSE57338/ | 4.2G | ✅ 已分析 |
| GSE116250 | 01_data/01_raw_geo/GSE116250/ | 20M | ✅ 已分析 |
| GSE55296 | 01_data/01_raw_geo/GSE55296/ | 1.8M | ✅ 已分析 |
| GDS4772 | 01_data/01_raw_geo/GDS4772/ | 7.0M | ⚠️ 分析阻塞(GPL超时) |
| GSE109816 | 01_data/01_raw_geo/GSE109816/ | 35M | ⏳ 后台运行中 |
| GSE183852 | 01_data/01_raw_geo/GSE183852/ | 17G | ✅ 已归档 |
| PXD010154 | 01_data/05_proteomics_pride/PXD010154/ | 3.1G | ✅ 已验证 |
| eQTLGen | 01_data/03_mr_data/eQTLGen/ | 309M | ✅ 已提取 |
| GTEx | 01_data/03_mr_data/gtex/ | 3.2G | ✅ 已提取 |
| HERMES v1 | 01_data/03_mr_data/gwas_summary/hermes_v1/ | 315M | ✅ 已归档 |
| HERMES v2 EUR | 01_data/03_mr_data/gwas_summary/hermes_v2_eur/ | 4.4G | ✅ 已归档 |
| HERMES v2 ALL | 01_data/03_mr_data/gwas_summary/hermes_v2_all/ | 3.2G | ✅ 已归档 |

## 今日关键教训

1. **伪重复是Visium最大陷阱**: 16004 spots ≠ 16004样本，有效n=5心脏
2. **健康异质性不可忽视**: GT_IZ_P13(64%零值) vs P15(10%零值)
3. **方法学警示 = 科学贡献**: "不显著"与"显著"同等重要
4. **并行策略有效**: 前台分析 + 后台下载，最大化时间产出
5. **及时止损**: GDS4772/Coloc阻塞时，立即转向替代方案

## 论文叙事最终版 (V61)

> "NDUFB7, a 137-aa microenvironment-sensitive accessory subunit of mitochondrial 
> Complex I, shows no robust overall downregulation in heart failure across four 
> independent platforms (Affymetrix n=330, RNA-seq n=101, p=0.52-0.81), but reveals 
> etiology-specific gradients (DCM>Ischemic, GSE55296 p=0.046) and substantial 
> inter-patient heterogeneity (10-64% zero-value spots across healthy hearts). 
> Spatial pseudoreplication correction attenuates the fibrotic-zone signal 
> (LMM: β=-0.148, p=0.81, n=5 hearts), positioning NDUFB7 as a hypothesis 
> requiring larger cohort validation rather than an established biomarker."

