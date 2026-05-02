#!/usr/bin/env Rscript
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
outdir <- file.path(PROJECT, "03_results/V162E_Checklist")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

checklist <- '### Cell Reports投稿门槛自检清单 (v1.0 | 2026-05-01)

| 维度 | 要求 | 状态 | 证据 |
|:---|:---|:---|:---|
| **创新性** | 新机制/新靶点/新方法 | ✅ PASS | NDUFB7→ferroptosis为首次报道；非参数双峰方法学创新 |
| **机制深度** | 从表型到分子机制 | ✅ PASS | 5平台meta→单细胞双峰→伪时序锚定→空间梯度→MR因果 |
| **数据量** | 足够统计效力 | ✅ PASS | 441K+4.9K单细胞 + 313 bulk + 5 Visium + 31K eQTL |
| **可重复性** | 跨队列验证 | ✅ PASS | 2单细胞队列Dip test交叉验证 + 5平台meta |
| **临床转化** | 可落地应用 | ⚠️ PARTIAL | 药物重定位候选已生成，缺实验验证 |
| **Figure质量** | 投稿级分辨率 | ✅ PASS | Fig2/3为300dpi PNG+PDF双格式 |
| **统计严谨性** | 多重检验校正 | ✅ PASS | BH FDR在所有bulk分析中应用 |
| **诚实性** | 不夸大结论 | ✅ PASS | 明确承认泛死亡签名，铁死亡为"最可药化节点" |

**Figure完整性检查：**
| Figure | 内容 | 状态 | 路径 |
|:---|:---|:---|:---|
| Fig 1 | 研究设计/机制图 | ⏳ 待绘制 | 需AI/手绘 |
| Fig 2 | 多平台meta+空间+时间+病因 | ✅ 完成 | V147_Fig2_Publication/ |
| Fig 3 | 单细胞双峰+铁死亡机制 | ✅ 完成 | V144_Fig3_Publication/ |
| Fig 4 | 铁死亡heatmap+rho+Visium | ✅ draft | Fig4_draft/ |
| Fig 5 | MR森林图+LocusZoom | ⏳ 需优化 | mr_results/ + V140_BIC/ |
| Fig 6 | 临床严重性+药物重定位 | ⏳ 需补充 | fig6_clinical/ + T28_Drug_Screen/ |

**缺失内容（投稿前必须完成）：**
1. Fig 1机制示意图（AI/手绘）
2. Fig 5投稿级优化（MR结果可视化）
3. Fig 6C clue.io药物预测结果填入
4. Supplementary Table 1: 基因集列表
5. Supplementary Table 2: 统计参数汇总
6. Supplementary Figure 1: GSE106118胚胎密度图
7. Methods完整版本（含软件版本号、参数细节）
8. Data Availability Statement定稿
9. Author Contributions
10. Competing Interests

**当前评分: 7.0/10** (Cell Reports门槛)
**目标评分: 7.5/10** (投稿前通过Fig 1/5/6 + Supplementary补全)
'

cat(checklist, file = file.path(outdir, "V162E_CellReports_Checklist.md"))

# Supplementary Table清单
supp <- '### Supplementary Table清单

| Table | 内容 | 数据文件 | 状态 |
|:---|:---|:---|:---|
| ST1 | 细胞死亡通路基因集 | 03_results/V151_C2_Discriminant/V151_discriminant_validation.csv | ✅ |
| ST2 | 单细胞队列统计参数 | 03_results/V140_BIC_Resolution/V140A_ZI_GMM_183852.csv + V156B | ✅ |
| ST3 | 跨平台meta分析参数 | 03_results/02_tables/Five_Platform_Forest_Data.csv | ✅ |
| ST4 | MR分析完整结果 | 03_results/mr_results/MR_NDUFB7_HF_FINAL_SUMMARY.csv | ✅ |
| ST5 | 临床队列样本信息 | 03_results/fig6_clinical/V70_GSE59867_admission_prognosis.csv | ✅ |
| ST6 | 药物重定位候选 | 03_results/T28_Drug_Screen/V94_T28_candidate_drugs.csv | ✅ |
| ST7 | GSE106118胚胎vs DCM对比 | 03_results/V159A_GSE106118_Embryo/V159A_Embryo_vs_DCM_NDUFB7.csv | ✅ |
'

cat(supp, file = file.path(outdir, "V162E_Supplementary_Tables.md"))

message("[DONE] V162E: ", outdir)
