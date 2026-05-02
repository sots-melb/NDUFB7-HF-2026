#!/bin/bash
cd "$(dirname "$0")"
PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"
SCRIPT_DIR="$PROJECT/02_scripts/01_geo_download"
LOG="$PROJECT/04_logs"
mkdir -p "$LOG"

echo ""
echo "╔══════════════════════════════════════════════════════════════════════════════╗"
echo "║     NDUFB7_Mito_2026 - 48小时高速GEO批量下载启动                             ║"
echo "╚══════════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "📊 磁盘空间:"
df -h "$PROJECT" | awk 'NR==2 {printf "   可用: %s (已用%s/%s)\n", $4, $3, $2}'
echo ""

echo "🚀 启动下载..."
cd "$SCRIPT_DIR"
nohup bash 05_final_download.sh > "$LOG/final_download.log" 2>&1 &
MAIN_PID=$!
echo $MAIN_PID > "$LOG/main.pid"
echo "   ✅ 主进程已启动 (PID: $MAIN_PID)"
echo ""
echo "📋 监控命令:"
echo "   实时监控: bash $SCRIPT_DIR/01_monitor.sh"
echo "   查看日志: tail -f $LOG/final_download.log"
echo "   数据校验: bash $SCRIPT_DIR/04_verify_data.sh"
echo ""
sleep 2
bash "$SCRIPT_DIR/01_monitor.sh"
