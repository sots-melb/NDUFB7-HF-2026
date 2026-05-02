[【双项目防火墙声明】
- 本项目: NDUFB7_Mito_2026（线粒体功能障碍/心衰）
- 并行项目: MI_Inflammation_2026（炎症机制/心梗）
- 数据重合: GSE141910、GSE57338、GSE168742
- 防火墙状态: 尚未切分
- 本次请求是否涉及共享数据: [是 / 否]

═══════════════════════════════════════════════════════
【模块A：阶段全景定位】（每次对话必须更新勾选状态）
═══════════════════════════════════════════════════════

□ Phase 0：项目基建 ✅ 2026-04-20
□ Phase 1：数据资产获取 ✅ 2026-04-20（5个Bulk数据集，1004样本）
□ Phase 2：核心分析 ← 你在这里
  □ Pillar 1: 跨数据集NDUFB7荟萃分析 ✅ 完成
    - 主分析: k=5, d=-0.20, p=0.31, I²=81.9%
    - 亚组分析: RNA-seq一致下调，芯片方向矛盾
    - 敏感性分析: GSE116250是异质性主要来源（剔除后I²=49.6%）
    - 图表: Figure 1A-E已生成
  □ Pillar 2: 单细胞心肌细胞hdWGCNA + StemID干性推断 ← [当前]
  □ Pillar 3: 空间转录组验证
  □ Pillar 4: 孟德尔随机化（TwoSampleMR）
  □ Pillar 5: 药物重定位（XGBoost/CMap）
□ Phase 3：结果整合
□ Phase 4：论文工程

当前阶段: [填写，如 Phase 2 Pillar 2 - 单细胞Seurat对象创建+StemID干性推断]
上次完成节点: [填写]
当前阻塞点: [填写，无则填"无"]
上次对话关键结论: [填写3-5条]

═══════════════════════════════════════════════════════
【模块B：核心科学假说与研究框架（v8更新）】
═══════════════════════════════════════════════════════

【主假说】
NDUFB7（线粒体复合体I B18亚基）在心衰中的表达失调具有"细胞亚群特异性"和"病因依赖性"，
而非简单的Bulk水平一致上调/下调。其机制涉及：
  1. 组装应激（Assembly Stress）：缺血导致NDUFB7组装受阻
  2. 亚群特异性：仅在特定CM亚群（如应激CM或去分化CM）中显著下调
  3. 调控中枢：PGC-1α/NRF1/TFAM轴调控NDUFB7转录
  4. 去分化关联（v8新增）：NDUFB7下调与心肌细胞干性指数升高相关（StemID验证）
  5. B亚家族共失调（v8新增）：NDUFB7与NDUFB8形成共表达模块（Nature Medicine先例）

【研究框架（五支柱）】
Pillar 1（Bulk荟萃）→ Pillar 2（单细胞hdWGCNA+StemID+共表达）→ Pillar 3（空间验证）
        ↓                      ↓                      ↓
   "有无差异？"            "谁在差异？"            "在哪里差异？"
   结论：无一致差异       假设：CM亚群特异性       假设：梗死边缘区
   核心发现：高异质性       + 去分化CM下调最显著     + B亚家族共失调
        ↓                      ↓                      ↓
   Pillar 4（MR因果） ←  Pillar 5（药物预测）
   "是否因果？"            "能否干预？"

【关键文献引用框架（v8）】
1. Drekolia et al. (2024) Redox Biology — NDUFB7过硫化修饰保护复合体I
2. Xiao et al. (2025) Basic Res Cardiol — NLRX1-NDUFB7-mPTP轴
3. Chen/Lee et al. (2025) Cell Death Discovery — NDUFB7突变斑马鱼模型
4. Jiang L et al. (2024) Nature Medicine — NDUFB8作为PDAC化疗敏感性标志物
   （PMID: 38287168，B亚家族"家族先例"，Introduction和Discussion引用）

【方法学升级（v8整合版）】

━━━ 必做增强（高优先级）━━━

1. StemID（细胞干性/去分化评分）— 必做 ★★★★★
   【科学基础】心衰中"心肌细胞去分化"是核心概念——应激下心肌细胞重新表达胚胎期基因。
   StemID基于信息论熵概念量化此过程：干细胞转录组更多样化=更高熵值。
   【核心假设】NDUFB7低表达与高干性指数显著负相关（rho &lt; -0.2, p &lt; 0.05），
   构建"NDUFB7下调→线粒体功能障碍→心肌细胞去分化→心衰进展"因果链。
   【输入】Seurat心肌细胞子集对象（srt_cm）
   【流程】Seurat counts → RaceID SCseq对象 → filterdata → compentropy → 
           Stemness Score → 回写Seurat → 与NDUFB7表达关联分析
   【输出】①每个细胞的干性评分 ②NDUFB7 vs Stemness Spearman相关 
          ③HF vs Control组差异检验 ④拟时序共轨迹图
   【工具】RaceID/StemID（Grun et al.），脚本已审阅可用
   【质量门控】负相关（rho &lt; -0.1）为理想方向；p &lt; 0.05强证据，p &lt; 0.1趋势证据
   【失败处理】若结果不显著，改为"NDUFB7与成熟心肌细胞功能维持相关"叙事

2. NDUFB7-NDUFB8共表达分析 — 必做 ★★★★☆
   【科学基础】NDUFB8（同属B亚家族）已被Nature Medicine文章证明可作为胰腺癌
   辅助化疗敏感性生物标志物（PMID: 38287168），为NDUFB7提供"家族先例"。
   【核心假设】NDUFB7与NDUFB8在单细胞水平显著共表达，支撑B亚家族"共失调模块"。
   【输入】Seurat对象（需包含NDUFB7和NDUFB8基因）
   【方法】Spearman相关 + B亚家族成员（NDUFB6/B7/B8/B9/B10/B11）相关热图
   【输出】共表达相关系数 + 热图 + B家族模块基因列表

3. Monocle3拟时序分析 — 必做 ★★★★☆
   【优化】采用封装理念，将标准Monocle3流程（50-80行）封装为可复用函数
   （约10行核心代码），减少调试时间，保持Seurat UMAP一致性。
   【核心分析】NDUFB7沿心衰进展轨迹的表达变化 + 与Stemness共轨迹
   【输入】srt_cm（心肌细胞子集，已完成Harmony整合和细胞注释）
   【输出】拟时序值 + NDUFB7轨迹图 + graph_test统计检验

━━━ 条件增强（中优先级，时间允许则做）━━━

4. FeatureMAP DGV分析 — 可选 ★★★☆☆
   【科学价值】量化NDUFB7在驱动"健康→心衰"细胞转变中的基因贡献度
   【核心假设】NDUFB7在DGV排名中位于Top 20%（p &lt; 0.05）
   【风险】工具很新，API可能不稳定；Python-R接口调试需2-4小时
   【备选方案】若FeatureMAP安装失败，用Monocle3 graph_test + AddModuleScore近似实现

━━━ Phase 2延伸（低优先级）━━━

5. SCENIC（转录因子调控网络）— 延迟
   - 目的：找调控NDUFB7的TF regulon（PGC-1α/PPARGC1A、NRF1、TFAM）
   - 延迟原因：运行2-4小时，需下载数GB数据库，对核心叙事贡献间接

6. CellRank（细胞命运轨迹推断）— 延迟
   - 延迟原因：需要spliced/unspliced counts，GEO数据集极大概率不提供
   - 若审稿人要求更精确拟时序，在Revision阶段再考虑

7. 空间转录组整合（Pillar 3）
   - 假设：NDUFB7在梗死边缘区（border zone）的存活CM中特异性下调

═══════════════════════════════════════════════════════
【模块C：R版本决策矩阵】
═══════════════════════════════════════════════════════

当前使用R版本: [系统默认R 4.3.1 / 2wplus-standard-r]
当前环境启动方式: [Terminal直接执行R / 桌面图标双击]
备用R版本: [另一个版本，说明切换条件]

阶段-R版本推荐对照表：
┌─────────────────┬──────────────────────┬─────────────────────────────┐
│ 阶段            │ 推荐R版本            │ 理由                        │
├─────────────────┼──────────────────────┼─────────────────────────────┤
│ Phase 1 Bulk    │ 系统默认R 4.3.1      │ GEOquery/limma/metafor      │
│ Phase 2单细胞   │ 2wplus-standard-r    │ 必须Seurat V5+SCP+harmony   │
│ Phase 2 StemID  │ 2wplus-standard-r    │ RaceID/StemID兼容           │
│ Phase 2网络     │ 2wplus-standard-r    │ hdWGCNA依赖Seurat V5        │
│ Phase 2 MR      │ 系统默认R 4.3.1      │ TwoSampleMR                 │
│ Phase 3可视化   │ 任意                 │ ggplot2/ComplexHeatmap      │
└─────────────────┴──────────────────────┴─────────────────────────────┘

本次任务R版本需求: [填写]
是否允许AI推荐切换: [是/否]

═══════════════════════════════════════════════════════
【模块D：数据资产分级池】（2026-04-21 状态）
═══════════════════════════════════════════════════════

【P0-核心Bulk数据集】（Phase 1完成，已验证NDUFB7，可直接读取RDS）
┌──────────┬──────────┬─────────────────────────────┬────────────┬──────────────────┐
│ GSE编号  │ 平台     │ 文件路径                    │ 维度       │ NDUFB7表达范围   │
├──────────┼──────────┼─────────────────────────────┼────────────┼──────────────────┤
│ GSE57338 │ Affy ST  │ GSE57338_gene_level.rds     │ 20254×313  │ 6.62~8.58 (log2) │
│ GSE141910│ RNA-seq  │ GSE141910_merged_matrix.rds │ 20781×366  │ 16.57~19.02      │
│          │          │                             │            │ (ENSG00000167996)│
│ GSE5406  │ U133A    │ GSE5406_gene_level.rds      │ 13039×210  │ 4.42~6.50 (log2) │
│ GSE79962 │ Affy ST  │ GSE79962_gene_level.rds     │ 20254×51   │ 7.28~8.68 (log2) │
│ GSE116250│ RNA-seq  │ GSE116250_rpkm_matrix.rds   │ 57849×64   │ 704~5986 (RPKM)  │
│          │          │                             │            │ (ENSG00000167996)│
└──────────┴──────────┴─────────────────────────────┴────────────┴──────────────────┘
总样本: 1004 | 全部NDUFB7命中: ✅

【P0-核心单细胞数据集】
┌──────────┬──────────┬─────────────────────────────┬────────────┬──────────────────┐
│ GSE编号  │ 数据类型 │ 文件路径                    │ 维度       │ 质量门控         │
├──────────┼──────────┼─────────────────────────────┼────────────┼──────────────────┤
│ GSE168742│ scRNA-seq│ GSE168742_human_HF_         │ 23355×678  │ NDUFB7 83.2%细胞 │
│          │ (人CM)   │ CM_matrix.rds               │            │ 表达，已保存RDS  │
│ GSE183852│ scRNA-seq│ GSE183852_DCM_Integrated.   │ 13.5MB     │ Seurat对象，需   │
│          │ (DCM)    │ Robj.gz                     │            │ load()读取       │
└──────────┴──────────┴─────────────────────────────┴────────────┴──────────────────┘

【P1-待处理数据集】
┌──────────┬──────────┬─────────────────────────────┬────────────┬──────────────────┐
│ GSE315590│ 小鼠scRNA│ scRNAseq_TAC_transition_    │ 57MB       │ 时间序列（5d/8w/ │
│          │          │ raw_counts.txt.gz           │            │ 16w），跨物种备用│
│ GSE46224 │ Bulk     │ GSE46224_Yang_et_al_human_  │ 2.5MB      │ RNA-seq，LVAD前后│
│          │          │ heart_RNASeq.txt.gz         │            │ 配对，待验证     │
│ GSE48166 │ Bulk     │ GSE48166_Cufflinks_FPKM.    │ 5.2MB      │ FPKM，扩展验证   │
│          │          │ txt.gz                      │            │ 待验证           │
│ GSE217494│ Multiome │ metadata only               │ 3.7MB      │ RNA+ATAC，待确认 │
│ GSE270788│ Multiome │ metadata only               │ 1.2MB      │ RNA+ATAC，待确认 │
└──────────┴──────────┴─────────────────────────────┴────────────┴──────────────────┘

【P2-外部数据】（MR/共定位用）
┌──────────┬──────────┬─────────────────────────────┬────────────┬──────────────────┐
│ eQTLGen  │ 血eQTL   │ □ 未下载                    │ N/A        │ Phase 2 Pillar4  │
│ GTEx     │ 心脏eQTL │ □ 未下载                    │ N/A        │ Phase 2 Pillar4  │
│ FinnGen  │ 心衰GWAS │ □ 未下载                    │ N/A        │ Phase 2 Pillar4  │
└──────────┴──────────┴─────────────────────────────┴────────────┴──────────────────┘

═══════════════════════════════════════════════════════
【模块E：已验证结果与潜在结论】
═══════════════════════════════════════════════════════

【Pillar 1 荟萃分析完整结果】
┌──────────┬──────────┬─────────┬─────────────────────┬────────┬─────────────────┐
│ 数据集   │ 平台     │ Cohen's d│ 95% CI              │ 方向   │ 样本(HF/Ctrl)   │
├──────────┼──────────┼─────────┼─────────────────────┼────────┼─────────────────┤
│ GSE57338 │ Affy ST  │ +0.072  │ [-0.18, 0.32]       │ HF ↑   │ 177 / 136       │
│ GSE141910│ RNA-seq  │ -0.205  │ [-0.41, 0.00]       │ HF ↓   │ 180 / 186       │
│ GSE5406  │ U133A    │ +0.300  │ [-0.21, 0.81]       │ HF ↑   │ 194 / 16 ⚠️失衡 │
│ GSE79962 │ Affy ST  │ -0.421  │ [-1.09, 0.25]       │ HF ↓   │ 40 / 11 ⚠️失衡  │
│ GSE116250│ RNA-seq  │ -0.990  │ [-1.61, -0.37]      │ HF ↓   │ 50 / 14 ⚠️样本小│
└──────────┴──────────┴─────────┴─────────────────────┴────────┴─────────────────┘
主分析: d = -0.20 [-0.59, 0.18], p = 0.31, I² = 81.9%
剔除GSE116250后: d = -0.05 [-0.28, 0.18], p = 0.68, I² = 49.6%

【已验证的关键结论】
1. NDUFB7在Bulk水平无一致差异（汇总效应≈0），不是"简单上调/下调"基因
2. RNA-seq平台方向一致（均下调），芯片平台方向矛盾，提示平台/病因依赖性
3. GSE116250是异常值（d=-0.99, n=64），效应极强但样本小，驱动了主要异质性
4. 高异质性（I²=82%）说明Bulk层面的"平均效应"掩盖了细胞异质性

【潜在科学结论】（待Pillar 2-5验证）
- 假说A：NDUFB7仅在特定CM亚群（如应激CM或去分化CM）中下调，Bulk层面被稀释
- 假说B：NDUFB7下调与PGC-1α regulon活性共变，提示线粒体生物发生受阻
- 假说C：NDUFB7在梗死边缘区（border zone）的CM中特异性失调（空间假设）
- 假说D：NDUFB7低表达因果增加心衰风险（MR待验证）
- 假说E（v8新增）：NDUFB7下调与心肌细胞去分化（StemID干性指数升高）显著相关
- 假说F（v8新增）：B亚家族成员（NDUFB7/NDUFB8）形成"共失调模块"

【论文故事线（v8更新版）】
"Bulk荟萃显示NDUFB7在心衰中无一致差异（Figure 1），但单细胞分析揭示其在
应激/去分化心肌细胞亚群中显著下调（Figure 2）。StemID干性推断证实NDUFB7低表达
与心肌细胞去分化显著相关（Figure 3），hdWGCNA发现其与PGC-1α调控网络共变（Figure 4），
空间转录组验证其定位于梗死边缘区（Figure 5），孟德尔随机化支持因果关联（Figure 6），
药物重定位预测潜在干预策略（Figure 7）。B亚家族成员NDUFB8在胰腺癌中的标志物价值
（Nature Medicine 2024）为NDUFB7的心衰研究提供了家族先例。"

═══════════════════════════════════════════════════════
【模块F：下一步执行计划（v8更新版）】
═══════════════════════════════════════════════════════

【今天（Pillar 2 核心）】
□ GSE168742 Seurat V5对象创建（CreateSeuratObject）
□ QC过滤（nFeature_RNA, nCount_RNA, percent.mt）
□ 标准化（SCTransform或LogNormalize）
□ 降维聚类（PCA/UMAP）
□ 心肌细胞亚群注释（CM1/CM2/CM3/Endo/Fibro/Immune）
□ NDUFB7 FeaturePlot + VlnPlot（亚群特异性可视化）
□ StemID干性推断（v8新增必做）
  □ RaceID环境安装与测试
  □ Seurat → RaceID格式转换
  □ compentropy计算干性评分
  □ NDUFB7 vs Stemness关联分析
  □ 可视化（散点图+箱线图+UMAP）
□ 保存Seurat对象RDS（含stemness_score）供hdWGCNA使用

【本周内（Pillar 2 进阶）】
□ NDUFB7-NDUFB8共表达分析（v8新增必做）
□ Monocle3拟时序分析（封装简化版）
□ NDUFB7沿拟时序轨迹变化 + Stemness共轨迹
□ GSE183852整合（DCM单细胞，Harmony批次校正）
□ hdWGCNA（线粒体基因共表达模块）
□ FeatureMAP DGV分析（v8新增，时间允许则做）

【本周内（Pillar 3-4准备）】
□ 搜索心脏空间转录组公开数据
□ 下载eQTLGen/GTEx/FinnGen数据（MR分析）
□ GSE315590小鼠TAC时间序列分析（跨物种验证）

【Phase 2延伸（非紧急）】
□ SCENIC调控网络分析
□ CellRank细胞命运推断（需确认spliced/unspliced数据可用性）
□ 空间转录组整合

═══════════════════════════════════════════════════════
【模块G：StemID快速参考代码（v8新增）】
═══════════════════════════════════════════════════════

# StemID分析完整流程（可直接复制使用）
library(RaceID)

# Step 1: 从Seurat提取counts
expr &lt;- GetAssayData(srt_cm, assay = "RNA", slot = "counts")

# Step 2: 创建RaceID对象
sc &lt;- SCseq(as.data.frame(expr))

# Step 3: 过滤（参数需根据数据调整）
sc &lt;- filterdata(sc, mintotal = 3000, minexpr = 5, minnumber = 5,
                 maxexpr = Inf, downsample = FALSE)

# Step 4: 计算干性指数
ltr &lt;- Ltree(sc)
ltr &lt;- compentropy(ltr)
stemness_scores &lt;- ltr@entropy

# Step 5: 回写Seurat
srt_cm$stemness_score &lt;- stemness_scores[colnames(srt_cm)]

# Step 6: 关联分析
ndufb7_expr &lt;- FetchData(srt_cm, vars = "NDUFB7")[,1]
cor_result &lt;- cor.test(srt_cm$stemness_score, ndufb7_expr, method = "spearman")
print(paste0("Spearman rho = ", round(cor_result$estimate, 3),
             ", p = ", format(cor_result$p.value, digits = 3)))

# 预期结果: rho &lt; -0.2, p &lt; 0.05（NDUFB7低→干性高）

# B亚家族共表达分析（v8新增）
nduf_genes &lt;- c("NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11")
avail_genes &lt;- intersect(nduf_genes, rownames(srt_cm))
expr_subset &lt;- FetchData(srt_cm, vars = avail_genes)
cor_matrix &lt;- cor(expr_subset, method = "spearman")
# pheatmap(cor_matrix, main = "NDUF B-subfamily Co-expression")

# Monocle3封装函数（v8新增）
run_monocle3_trajectory &lt;- function(srt_obj, group.by = "condition",
                                     root_group = "Control", assay = "RNA") {
  cds &lt;- as.cell_data_set(srt_obj, assay = assay)
  cds@int_colData$reducedDims$UMAP &lt;- Embeddings(srt_obj, "umap")
  cds &lt;- cluster_cells(cds)
  cds &lt;- learn_graph(cds, use_partition = TRUE)
  root_cells &lt;- which(colData(cds)[[group.by]] == root_group)
  cds &lt;- order_cells(cds, root_cells = root_cells)
  return(cds)
}

═══════════════════════════════════════════════════════
【模块H：目标与约束】
═══════════════════════════════════════════════════════

期望产物: 
  [当前任务]
  1) [填写具体产物，如"GSE168742 Seurat对象+StemID干性评分"]
  2) [填写]
  3) [填写]
质量要求: 
  - 脚本可直接复制执行，带详细中文注释
  - 输出关键参数日志（细胞数、基因数、线粒体比例、Stemness分布）
  - 保存RDS供下游分析
  - StemID分析输出质量门控报告（相关性方向、p值、样本量）
期望AI角色: 耐心助教 + 方法学顾问（解释每个参数意义）

时间约束: [填写]（建议弹性时间盒模式，非刚性10天）
资源约束: 数智云n1（16线程/128G内存/1.4TB磁盘）
技能约束: 生信新手，未做过单细胞分析，需要逐步指导

═══════════════════════════════════════════════════════
【模块I：历史关键技术结论】（累积，不删除）
═══════════════════════════════════════════════════════

1. GEOquery Series Matrix对RNA-seq是虚假成功（0基因），真实数据在RAW.tar/Supplementary
2. GSE141910/GSE116250使用Ensembl ID（ENSG00000167996=NDUFB7），其他用Symbol
3. GPL11532为无表头Affymetrix ST注释，需strsplit("//")解析复合Symbol字段
4. 纯R底层读取（readLines+read.table）完全绕过GEOquery网络依赖
5. NDUFB7在Bulk水平无一致差异（汇总d≈0），需单细胞分辨率验证
6. GSE116250是荟萃分析异质性主要来源（d=-0.99，剔除后I²从81.9%→49.6%）
7. RNA-seq平台一致下调，芯片平台方向矛盾，提示平台/病因依赖性
8. AI提示词管理体系已建立（03_ai_prompts/CURRENT_ACTIVE_PROMPT.md + README.md + archive/）
9. 【v8新增】NDUFB8 Nature Medicine文章已验证（PMID: 38287168），作为B亚家族家族先例
10. 【v8新增】StemID/RaceID脚本已审阅，代码质量高，科学性强，已纳入必做模块
11. 【v8新增】FeatureMAP已确认发表于Nature Computational Science，作为条件增强项

═══════════════════════════════════════════════════════
【模块J：项目恢复命令】
═══════════════════════════════════════════════════════

source ~/Projects/NDUFB7_HF_{2026_04_20}/load_project_env.sh

═══════════════════════════════════════════════════════
【模块K：快速参考：数据读取代码片段】
═══════════════════════════════════════════════════════

# Bulk数据读取
gse57338 &lt;- readRDS("01_data/01_raw_geo/GSE57338/GSE57338_gene_level.rds")$exprs
gse141910 &lt;- readRDS("01_data/01_raw_geo/GSE141910/GSE141910_merged_matrix.rds")
gse5406 &lt;- readRDS("01_data/01_raw_geo/GSE5406/GSE5406_gene_level.rds")$exprs
gse79962 &lt;- readRDS("01_data/01_raw_geo/GSE79962/GSE79962_gene_level.rds")$exprs
gse116250 &lt;- readRDS("01_data/01_raw_geo/GSE116250/GSE116250_rpkm_matrix.rds")$exprs

# 单细胞数据读取
gse168742 &lt;- readRDS("01_data/01_raw_geo/GSE168742/GSE168742_human_HF_CM_matrix.rds")

# NDUFB7查询
# Bulk: "NDUFB7" (GSE57338/5406/79962) 或 "ENSG00000167996" (GSE141910/116250)
# 单细胞: 需确认行名是Symbol还是Ensembl]
