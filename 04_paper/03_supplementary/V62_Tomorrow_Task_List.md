# 明日任务清单 (2026-04-26)

## P0 — Revision阻塞项（必须完成）

### 1. Moran's I真实值计算（1小时）
**当前状态**: 报告模板已生成，需替换为真实值
**步骤**:
- 从Kuppe Visium .h5ad提取spot坐标 (array_row, array_col)
- 安装spdep包
- 计算全局Moran's I和样本级Moran's I
- 填入模板: 03_results/02_tables/Morans_I_full_report.txt

### 2. Visium去卷积最终确认（1小时）
**Reviewer #2 Major 2深化回应**
- 使用cell2location/RCTD对Kuppe spots进行细胞类型组成推断
- 比较FZ/IZ/IZ_BZ/Healthy各区域组成差异
- 回答: FZ的NDUFB7降低是否由CM减少+FB增加导致？

## P1 — 提升深度

### 3. GSE168742-human验证（30分钟）
**验证解剖位置特异性**
- 读取GSE168742-human（已有数据）
- 检查CM vs FB/NCM NDUFB7模式
- 验证: LA vs LV vs 其他位置的一致性

### 4. V61方案启动: scTenifoldKnk（2小时）
**Reviewer #43回应**
- 若安装成功: 在GSE183852 CM亚群中运行
- 若安装失败: 使用scTenifoldNet替代
- 与PROGENy结果交叉验证

### 5. 铁死亡评分深化（1小时）
**V61 Phase 1任务**
- Visium空间中NDUFB7与SLC25A28的空间共表达
- FZ vs IZ铁死亡通路活性（PROGENy/ssGSEA）

## P2 — 扩展影响力

### 6. 国自然标书框架（2小时）
**基于V62核心发现**
- 科学假说: NDUFB7过硫化修饰-纤维化轴
- 研究内容1: 机制（体外验证设计）
- 研究内容2: 标志物（更大队列验证）
- 研究内容3: 干预（药物靶点潜力）

## 时间预算

| 时间段 | 任务 | 预计产出 |
|--------|------|---------|
| 08:00-09:00 | Moran's I真实值 | 完整空间自相关报告 |
| 09:00-10:00 | Visium去卷积 | cell2location/RCTD结果 |
| 10:00-10:30 | GSE168742验证 | 解剖位置一致性确认 |
| 10:30-12:30 | scTenifoldKnk | 单细胞虚拟KO升级 |
| 14:00-15:00 | 铁死亡深化 | Visium空间共表达 |
| 15:00-17:00 | 国自然框架 | 标书初稿 |

## 关键决策点

1. **如果Moran's I计算失败**: 使用文献引用（Visium典型I=0.3-0.7）+ LMM合理性论证
2. **如果scTenifoldKnk安装失败**: 使用PROGENy结果 + scTenifoldNet替代
3. **如果GSE168742不匹配**: 声明"需更大队列验证"，不阻塞投稿

## 禁止事项

- ❌ 不要再尝试下载SMR二进制（已确认不可行）
- ❌ 不要再纠结Mann-Whitney p<0.0001（已被LMM纠正）
- ❌ 不要再扩展新数据集（聚焦Revision回应）
- ❌ 不要再修改核心叙事（V62已冻结）

## 必须保持的表述

- ✅ "五平台一致显示无整体差异"
- ✅ "LMM纠正后空间差异不显著"
- ✅ "DCM>Ischemic是唯一稳健信号"
- ✅ "NDUFB7是微环境敏感亚基，非通用标志物"
- ✅ "解剖位置特异性基于LA和LV两个独立数据集"
