#!/bin/bash
PROJECT_ROOT="$HOME/Projects/NDUFB7_HF_{2026_04_20}"
BACKUP_DIR="$PROJECT_ROOT/manuscript/supplementary_scripts"
mkdir -p "$BACKUP_DIR"

SP=$1
shift

if [ ! -f "$SP" ]; then
  echo "[错误] 脚本不存在: $SP"
  exit 1
fi

SN=$(basename "$SP")
TS=$(date +%Y%m%d_%H%M%S)
AF="$BACKUP_DIR/${TS}_${SN}"

{
  echo "# =============================================================================="
  echo "# 论文审计信息（自动生成）"
  echo "# 原始脚本: $SP"
  echo "# 审计时间: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "# 执行用户: $(whoami)"
  echo "# 服务器: $(hostname)"
  echo "# 工作目录: $(pwd)"
  echo "# R版本: $(R --version | head -1 2>/dev/null || echo N/A)"
  echo "# =============================================================================="
  cat "$SP"
} > "$AF"

echo "[审计] 已归档: $AF"

if [[ "$SP" == *.R ]]; then
  RUNNER="Rscript"
elif [[ "$SP" == *.sh ]]; then
  RUNNER="bash"
else
  echo "[错误] 仅支持.R和.sh"
  exit 1
fi

echo "[执行] $RUNNER $SP $@"
START=$(date +%s)
$RUNNER "$SP" "$@"
EC=$?
END=$(date +%s)
DUR=$((END - START))

echo "[执行] 耗时: ${DUR}秒 | 退出码: $EC"
exit $EC
