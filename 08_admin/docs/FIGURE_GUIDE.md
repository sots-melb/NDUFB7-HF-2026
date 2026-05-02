# NDUFB7_Mito_2026 论文Figure草图指南

## Figure 1: Bulk荟萃分析（已完成 ✅）
- **1A**: 森林图（k=5, d=-0.20, p=0.31, I²=81.9%）
- **1B**: 亚组分析（RNA-seq vs Microarray）
- **1C**: 敏感性分析（剔除GSE116250后I²=49.6%）
- **1D**: 漏斗图
- **1E**: 发表偏倚检验
- **文件**: 03_results/08_meta_analysis/ 目录下PDF
- **状态**: ✅ 已完成，可直接放入论文

## Figure 2: 单细胞NDUFB7亚群特异性（已完成 ✅）
- **2A**: UMAP分群（6 clusters, 503 cells, HF 419 + Control 84）
- **2B**: NDUFB7 FeaturePlot（显示Cluster 4极低表达）
- **2C**: NDUFB7 Violin by Cluster（Cluster 4: 0.36 vs Cluster 3: 3.00）
- **2D**: HF vs Control细胞比例堆叠图
- **2E**: Cluster 4特征基因热图
- **文件**: 
  - 03_results/09_single_cell/06_umap_condition.pdf
  - 03_results/09_single_cell/07_ndufb7_featureplot.pdf
  - 03_results/09_single_cell/08_ndufb7_violin.pdf
- **状态**: ✅ 已完成

## Figure 3: B亚家族共塌陷与去分化否定（已完成 ✅）
- **3A**: B亚家族（NDUFB6/7/8/9/10/11）Cluster热图
- **3B**: NDUFB7-NDUFB8散点图（rho=0.174, p=8.55e-05）
- **3C**: B亚家族相关性热图
- **3D**: StemID干性推断散点图（rho=+0.08, p=0.073 → 否定去分化）
- **3E**: Monocle3拟时序图（NDUFB7无显著趋势）
- **文件**:
  - 03_results/09_single_cell/14_B_family_cor_heatmap.pdf
  - 03_results/09_single_cell/15_NDUFB7_vs_NDUFB8_scatter.pdf
  - 03_results/09_single_cell/16_B_family_cluster_heatmap.pdf
  - 03_results/09_single_cell/09_NDUFB7_stemness_scatter.pdf
  - 03_results/09_single_cell/17_NDUFB7_pseudotime.pdf
- **状态**: ✅ 已完成

## Figure 4: OXPHOS模块塌陷（已完成 ✅，hdWGCNA为升级项）
- **4A**: 64基因OXPHOS模块相关性热图
- **4B**: 模块评分UMAP（Cluster 4为蓝色低分）
- **4C**: 模块评分Violin（Cluster 4 vs 其他HF: p=2.15e-17）
- **4D**: 模块基因Cluster热图（显示跨复合体协调下调）
- **4E**: HF vs Control模块评分箱线图（p=9.14e-31）
- **文件**:
  - 03_results/09_single_cell/30_ndufb7_module_heatmap.pdf
  - 03_results/09_single_cell/31_module_score_umap.pdf
  - 03_results/09_single_cell/32_module_score_violin.pdf
  - 03_results/09_single_cell/33_module_genes_cluster_heatmap.pdf
- **状态**: ✅ 已完成，简化模块足够发表
- **升级**: ⏳ hdWGCNA hub基因网络图（若17:57后安装成功则加入4F）

## Figure 5: 空间转录组验证（待执行 ⏳）
- **5A**: 心脏组织H&E + Visium spot overlay
- **5B**: NDUFB7空间表达热图（预期梗死边缘区下调）
- **5C**: 病理分区比较（infarct vs border vs remote）
- **5D**: 空间模块评分（验证Cluster 4型细胞定位）
- **数据**: 待从GSE225295/GSE210867获取
- **状态**: ⏳ 候选清单已生成，待下载数据

## Figure 6: 孟德尔随机化（待执行 ⏳）
- **6A**: SNP效应散点图（NDUFB7 eQTL vs HF GWAS）
- **6B**: 森林图（各SNP Wald Ratio）
- **6C**: 漏斗图
- **6D**: 留一法敏感性分析
- **数据**: GTEx eQTL下载中（后台），FinnGen待下载
- **状态**: ⏳ 本地框架已保存，数据获取中

## Figure 7: 药物重定位（待执行 ⏳）
- **7A**: NDUFB7模块基因query signature
- **7B**: CMap/L1000候选药物热图
- **7C**: 药物-靶点网络图
- **状态**: ⏳ Pillar 5待启动

---

## 当前论文完成度评估
- Figure 1: 100% ✅
- Figure 2: 100% ✅
- Figure 3: 100% ✅
- Figure 4: 90% ✅（简化模块完成，hdWGCNA为+10%升级）
- Figure 5: 0% ⏳（数据获取阶段）
- Figure 6: 10% ⏳（框架就绪，数据下载中）
- Figure 7: 0% ⏳（未启动）

**结论**: Figure 1-4已足够支撑一篇完整的单细胞+荟萃分析论文（短篇/letter级别）。
Figure 5-7为扩展验证，完成后可冲击更高水平期刊。
