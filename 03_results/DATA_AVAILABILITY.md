# Data Availability Statement

## 已冻结核心数据（Tier S）

| 数据集 | 平台 | 样本量 | 状态 | 路径 |
|:-------|:-----|:-------|:-----|:-----|
| GSE183852 | snRNA-seq | 191,795 cells | ✅ 已分析 | 01_data/02_processed/GSE183852/ |
| Kuppe Visium | Visium V1 | 6 hearts, 19,525 spots | ✅ 已分析 | 01_data/02_spatial/ |
| GSE57338 | Affymetrix 1.1 ST | 313 LV samples | ✅ 已分析 | 01_data/01_raw_geo/GSE57338/ |
| GSE116250 | RNA-seq RPKM | 64 samples | ✅ 已分析 | 01_data/01_raw_geo/GSE116250/ |
| GSE55296 | RNA-seq count | 37 samples | ✅ 已分析 | 01_data/01_raw_geo/GSE55296/ |
| PXD010154 | MaxQuant蛋白组 | 19 fractions | ✅ 已提取 | 05_proteomics_pride/ |

## 已下载待分析数据（Tier A）

| 数据集 | 平台 | 样本量 | 状态 | 路径 |
|:-------|:-----|:-------|:-----|:-----|
| GSE154170 | RNA-seq | 16 ICM samples | ⏳ 分组待修复 | 01_data/01_raw_geo/GSE154170/ |
| GSE59867 | Affymetrix | 28 HF samples | ⏳ 预后分析 | 01_data/01_raw_geo/GSE59867/ |
| GSE168742 | scRNA-seq | 637 cells | ⏳ 验证待执行 | 01_data/01_raw_geo/GSE168742/ |
| GSE315590 | scRNA-seq | 10,164 cells (mouse TAC) | ✅ 已分析 | 01_data/01_raw_geo/GSE315590/ |

## V79新下载数据（series matrix）

| 数据集 | 类别 | 状态 | 路径 |
|:-------|:-----|:-----|:-----|
| GSE226314, GSE185100, GSE288222, GSE302337, GSE214611 | Tier S | ✅ series matrix | 01_data/01_raw_geo/V79_20260430/Tier_S/ |
| GSE269705, GSE247468, GSE222144, GSE214731, GSE138262, GSE227734 | Tier A | ✅ series matrix | 01_data/01_raw_geo/V79_20260430/Tier_A/ |
| GSE314910, GSE310386, GSE316390, GSE317680, GSE300585, GSE271676, GSE290577 | Spatial | ✅ series matrix | 01_data/01_raw_geo/V79_20260430/Spatial/ |

## 数据筛选流程

系统检索75个GEO数据集 + 2个蛋白数据库 → 按以下标准筛选：
1. 人类心脏组织
2. NDUFB7可检测
3. 样本量>10
4. 包含病因/空间信息

最终纳入：7个核心数据集（4转录组 + 1空间 + 2蛋白验证）
筛选流程图见Supplementary Figure 1。

## 数据获取限制

部分数据集因以下原因未纳入：
- 物种非人类（排除）
- 样本量<10（排除）
- NDUFB7检测率<5%（排除）
- 数据格式不兼容且无法修复（排除）
- 原始数据未公开（记录为限制）

## 代码可用性

所有分析脚本已归档至GitHub（sots-melb/NDUFB7-HF-2026），包含：
- 数据处理脚本（R/Python）
- 统计分析代码
- Figure生成脚本
- 错误日志（E1-E110）

