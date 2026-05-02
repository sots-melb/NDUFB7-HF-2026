#!/bin/bash
PROJECT_ROOT="$HOME/Projects/NDUFB7_HF_{2026_04_20}"
OUT="$PROJECT_ROOT/04_docs/02_meetings/ai_prompt_auto_$(date +%Y%m%d_%H%M%S).md"

echo "=== 正在扫描项目状态 ==="

{
  echo "# NDUFB7项目AI提示词（自动生成 $(date '+%Y-%m-%d %H:%M:%S')）"
  echo ""
  echo "## 1. R版本"
  R --version | head -1 | sed 's/^/  /'
  echo ""
  echo "## 2. 数据资产状态（自动扫描）"
  for d in GSE141910 GSE57338 GSE116250 GSE5406 GSE46224 GSE48166 GSE55296 GSE79962 GSE168742; do
    dp="$PROJECT_ROOT/01_data/01_raw_geo/$d"
    if [ -d "$dp" ]; then
      sz=$(du -sh "$dp" 2>/dev/null | cut -f1)
      fc=$(ls -1 "$dp" 2>/dev/null | wc -l)
      echo "  * $d: $sz ($fc files)"
    else
      echo "  * $d: ❌ 目录不存在"
    fi
  done
  echo ""
  echo "## 3. 最近24小时修改的脚本"
  find "$PROJECT_ROOT/02_scripts" -type f \( -name "*.R" -o -name "*.sh" \) -mmin -1440 2>/dev/null | head -5 | sed 's/^/  /'
  echo ""
  echo "## 4. 最近日志尾部"
  LL=$(ls -t "$PROJECT_ROOT/05_logs"/*.log 2>/dev/null | head -1)
  if [ -n "$LL" ]; then tail -5 "$LL" | sed 's/^/  /'; fi
  echo ""
  echo "## 5. 预填充AI提示词模板（请复制第5节到AI窗口）"
  echo ""
  echo '【双项目防火墙声明】'
  echo '- 本项目: NDUFB7_Mito_2026'
  echo '- 并行项目: MI_Inflammation_2026'
  echo '- 防火墙状态: 尚未切分'
  echo '- 本次请求是否涉及共享数据: 是'
  echo ''
  echo '【项目上下文】'
  echo '- 当前阶段: Phase 1 - 数据准备（请根据上方状态精确填写）'
  echo '- 当前任务: [待填写，如：读取GSE141910_RAW.tar并合并表达矩阵]'
  echo '- 上次进度: [见上方第4节日志]'
  echo '- 上次对话关键结论: GEOquery Series Matrix对RNA-seq是虚假成功'
  echo ''
  echo '【环境信息】'
  echo '- 服务器: 数智云独享计算型-n1 (16线程/128G内存/1.5T存储)'
  echo "- R版本: $(R --version | head -1)"
  echo '- 工作目录: ~/Projects/NDUFB7_HF_{2026_04_20}'
  echo '- 恢复命令: source ~/Projects/NDUFB7_HF_{2026_04_20}/load_project_env.sh'
  echo ''
  echo '【数据状态】'
  echo '[见上方第2节数据资产状态]'
  echo ''
  echo '【目标输出】'
  echo '- 期望产物: [待填写]'
  echo '- 质量要求: [待填写]'
  echo '- 期望AI角色: 耐心助教（生信新手）'
  echo ''
  echo '【限制条件】'
  echo '- 时间约束: [待填写]'
  echo '- 资源约束: 内存128G / 磁盘1.4TB'
  echo '- 技能约束: 生信新手，需要详细注释'
} > "$OUT"

echo ""
echo "╔════════════════════════════════════════════════════════════╗"
echo "║  AI提示词已自动生成                                        ║"
echo "╠════════════════════════════════════════════════════════════╣"
echo "║  文件: $OUT"
echo "║  查看: cat $OUT"
echo "║  复制第5节（预填充模板）粘贴到新AI窗口                     ║"
echo "╚════════════════════════════════════════════════════════════╝"
