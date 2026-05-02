# 明日任务清单 (2026-04-26)

## 🔴 P0 — 必须完成 (Revision阻塞项)

### 1. GDS4772补全 (30分钟)
**阻塞原因**: GEOquery GPL6244下载超时
**解决策略**: 
- 手动下载GPL6244注释文件: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6244
- 或直接使用GDS4772.soft中的探针ID暴力搜索NDUFB7
- 提取统计值后更新五平台森林图（替换占位符）

**验证标准**: 
- [ ] GDS4772_NDUFB7_stats.csv存在
- [ ] 森林图更新为真实五平台（非占位符）
- [ ] Reviewer #1 Major 2完全回应

### 2. GitHub代码上传 (1小时)
**Reviewer #3 Minor要求**
**步骤**:
1. 初始化GitHub仓库: `git init` / `git remote add origin`
2. 整理脚本到 `02_scripts/` 标准结构
3. 生成README.md（含依赖包版本）
4. 生成sessionInfo.txt（R包版本记录）
5. 推送至GitHub获取DOI（通过Zenodo集成）

**文件清单**:
- [ ] 所有R分析脚本（去卷积、混合效应模型、WGCNA、PROGENy）
- [ ] Python辅助脚本（h5ad提取、格式转换）
- [ ] bash下载/归档脚本
- [ ] 数据可用性声明
- [ ] 运行说明文档

## 🟡 P1 — 提升接收概率

### 3. GTEx eQTL正确文件获取 (30分钟)
**阻塞原因**: eQTLGen ENSG映射错误（映射到FTH1）
**解决策略**:
- 下载GTEx v8 Heart-LV正确eQTL文件
- 确认ENSG00000167996对应NDUFB7
- 重新运行Coloc共定位分析

### 4. GSE109816单细胞完成 (30分钟)
**当前状态**: 后台可能已完成
**检查**: `ls -lh 03_results/02_tables/GSE109816_NDUFB7_raw.rds`
**若完成**: 生成CM vs FB NDUFB7差异图

### 5. 论文Methods补全 (1小时)
**当前**: 已有草稿 (`03_results/paper_methods.md`)
**补全**:
- [ ] 混合效应模型详细参数（lme4公式、随机效应结构）
- [ ] 四平台统计方法统一描述
- [ ] 去卷积方法说明（cell2location引用）
- [ ] 遗传学分析流程（TwoSampleMR/Coloc）

## 🟢 P2 — 长期准备

### 6. 国自然标书框架 (2小时)
**基于V61核心发现**:
- 科学假说: NDUFB7过硫化修饰-纤维化轴
- 研究内容1: 机制（体外验证设计）
- 研究内容2: 标志物（更大队列验证）
- 研究内容3: 干预（药物靶点潜力）

### 7. 体外验证实验设计 (1小时)
**H9C2/AC16细胞系**:
- siNDUFB7 knockdown
- 检测: OXPHOS活性、ROS水平、纤维化标记物
- 对照: scramble siRNA

## 时间预算

| 时间段 | 任务 | 预计产出 |
|--------|------|---------|
| 08:00-08:30 | GDS4772补全 | 五平台森林图最终版 |
| 08:30-09:30 | GitHub上传 | 代码仓库+DOI |
| 09:30-10:00 | GTEx eQTL | Coloc重跑 |
| 10:00-10:30 | GSE109816 | 单细胞验证图 |
| 10:30-11:30 | Methods补全 | 完整Methods段落 |
| 14:00-16:00 | 国自然框架 | 标书初稿 |

## 关键决策点

1. **如果GDS4772仍失败**: 使用已有四平台森林图，在Figure legend中标注"GDS4772 pending"，不阻塞投稿
2. **如果GitHub上传失败**: 使用Zenodo直接上传tar.gz获取DOI
3. **如果GTEx eQTL仍有问题**: 使用eQTLGen全血作为exposure（跨组织MR，在Discussion中说明局限性）

## 禁止事项

- ❌ 不要再尝试下载SMR二进制（已确认不可行）
- ❌ 不要再纠结Mann-Whitney p<0.0001（已被LMM纠正）
- ❌ 不要再扩展新数据集（聚焦Revision回应）
- ❌ 不要再修改核心叙事（V61已冻结）

## 必须保持的表述

- ✅ "四平台一致显示无整体差异"
- ✅ "LMM纠正后空间差异不显著"
- ✅ "DCM>Ischemic是唯一稳健信号"
- ✅ "NDUFB7是微环境敏感亚基，非通用标志物"

