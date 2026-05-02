# AI提示词管理规范 v1.0

## 目录结构
03_ai_prompts/
├── CURRENT_ACTIVE_PROMPT.md    ← 当前活跃提示词（打开新窗口时复制此文件）
├── README.md                   ← 本文件
├── archive/                    ← 历史版本归档
│   ├── v4_20260420_phase1.md
│   ├── v5_20260420_phase2_pillar1.md
│   └── v6_20260421_phase2_pillar2.md
└── snippets/                   ← 常用片段（如R版本矩阵、数据资产表）
    ├── r_version_matrix.md
    └── data_assets_table.md

## 使用流程（每次新窗口）

1. **复制**：打开 `CURRENT_ACTIVE_PROMPT.md`，全选复制
2. **粘贴**：粘贴到新AI窗口的第一条消息
3. **更新**：根据当前任务填写【模块A】中的 `[填写]` 占位符
4. **对话**：正常进行技术分析
5. **归档**：如果对话中有重大更新（如新数据集验证、新方法突破），保存为新版本

## 版本命名规则
ai_prompt_v{主版本}_{YYYYMMDD}_{阶段标识}.md

主版本对应：
- v1-v9: Phase 1 数据准备
- v10-v19: Phase 2 核心分析  
- v20-v29: Phase 3 结果整合
- v30+: Phase 4 论文工程

## 更新触发条件（满足任一即更新）
□ 完成了新的里程碑（如Phase切换、Pillar完成）
□ 发现了关键技术结论（如"GSE116250是异质性来源"）
□ 数据资产发生重大变化（新增/删除数据集）
□ R版本切换或环境变更

## 禁止事项
✗ 不要在对话中直接修改CURRENT_ACTIVE_PROMPT.md（容易截断）
✗ 不要删除archive/中的历史版本（审计追踪）
✗ 不要把提示词散落在04_docs其他子目录
