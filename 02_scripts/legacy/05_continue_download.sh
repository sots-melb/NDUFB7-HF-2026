#!/bin/bash
PROJECT_ROOT="$HOME/Projects/NDUFB7_HF_{2026_04_20}"
DATA_RAW="$PROJECT_ROOT/01_data/01_raw_geo"
LOG_DIR="$PROJECT_ROOT/05_logs"
mkdir -p "$DATA_RAW" "$LOG_DIR"

LOG="$LOG_DIR/continue_download_$(date +%Y%m%d_%H%M%S).log"
echo "[$(date '+%H:%M:%S')] 继续下载缺失GSE（HTTPS模式）..." | tee -a "$LOG"

get_dir() {
    local num=${1#GSE}
    echo "GSE${num:0:3}nnn"
}

dl_one() {
    local gse=$1
    local target="$DATA_RAW/$gse"
    mkdir -p "$target"
    local sm="$target/${gse}_series_matrix.txt.gz"
    
    # 跳过已存在且完好的文件（>100KB）
    if [ -f "$sm" ] && gunzip -t "$sm" 2>/dev/null && [ $(stat -c%s "$sm") -gt 100000 ]; then
        echo "[$gse] ✅ 已存在 ($(du -sh "$sm" | cut -f1))" | tee -a "$LOG"
        return 0
    fi
    
    rm -f "$sm" "${sm}.tmp"
    local gdir=$(get_dir "$gse")
    local url="https://ftp.ncbi.nlm.nih.gov/geo/series/$gdir/$gse/matrix/${gse}_series_matrix.txt.gz"
    
    echo "[$gse] 下载..." | tee -a "$LOG"
    wget -c -t 3 --timeout=120 --progress=dot:giga -O "$sm.tmp" "$url" 2>>"$LOG" && \
    mv "$sm.tmp" "$sm" && \
    echo "[$gse] ✅ 完成 ($(du -sh "$sm" | cut -f1))" | tee -a "$LOG" || \
    echo "[$gse] ❌ 失败" | tee -a "$LOG"
    sleep 2
}

# 缺失GSE列表（根据之前资产确认）
for gse in GSE5406 GSE116250 GSE46224 GSE48166 GSE55296 GSE271946 GSE315590 GSE217494 GSE270788; do
    dl_one "$gse"
done

echo "[$(date '+%H:%M:%S')] 结束" | tee -a "$LOG"
