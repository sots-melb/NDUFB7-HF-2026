#!/bin/bash
REFRESH=${1:-10}
PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"
GEO="$PROJECT/01_data/01_raw_geo"
LOG="$PROJECT/04_logs"
G='\033[32m'; R='\033[31m'; Y='\033[33m'; C='\033[36m'; NC='\033[0m'; B='\033[1m'
trap 'echo -e "\n${G}监控已退出${NC}"; exit 0' INT
while true; do
    clear
    echo -e "${C}${B}╔══════════════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${C}${B}║          NDUFB7_Mito_2026 - 实时下载监控面板                                ║${NC}"
    echo -e "${C}${B}║          $(date '+%Y-%m-%d %H:%M:%S') | 刷新: ${REFRESH}s                              ║${NC}"
    echo -e "${C}${B}╚══════════════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    local avail_gb=$(df -BG "$PROJECT" | awk 'NR==2 {gsub(/G/,"",$4); print $4}')
    local used_pct=$(df -h "$PROJECT" | awk 'NR==2 {print $5}')
    local disk_color="$G"; [ "$avail_gb" -lt 5 ] && disk_color="$Y"; [ "$avail_gb" -lt 2 ] && disk_color="$R"
    echo -e "  ${B}磁盘空间:${NC} ${disk_color}${avail_gb}GB 可用${NC} (已用: $used_pct)"
    echo ""
    echo -e "  ${B}活跃下载进程:${NC}"
    local active_pids=$(pgrep -f "wget.*GSE\|wget.*GTEx\|wget.*FinnGen" 2>/dev/null)
    if [ -n "$active_pids" ]; then
        for pid in $active_pids; do
            local cmd=$(ps -p $pid -o args= 2>/dev/null | head -c 80)
            local cpu=$(ps -p $pid -o %cpu= 2>/dev/null | xargs)
            echo -e "    ${G}●${NC} PID:$pid | CPU:${cpu}% | ${cmd}..."
        done
    else
        echo -e "    ${Y}○ 无活跃下载进程${NC}"
    fi
    echo ""
    echo -e "  ${B}已下载数据:${NC}"
    for gse_dir in $(ls -1 "$GEO" 2>/dev/null | sort); do
        local file_count=$(ls -1 "$GEO/$gse_dir"/*.gz 2>/dev/null | wc -l)
        [ "$file_count" -gt 0 ] && printf "    📁 %-15s (%d files)\n" "$gse_dir" "$file_count"
    done
    echo ""
    echo -e "  ${C}────────────────────────────────────────────────────────────────────────${NC}"
    echo -e "  ${Y}Ctrl+C 退出 | 日志: tail -f $LOG/*.log${NC}"
    echo -e "  ${C}────────────────────────────────────────────────────────────────────────${NC}"
    sleep $REFRESH
done
