#!/bin/bash
# NDUFB7项目批量GEO下载脚本 v2.0
# 功能：断点续传、自动重试、日志审计、项目目录适配
# 用法：bash 00_batch_geo_download.sh

set -euo pipefail

PROJECT_ROOT="$HOME/Projects/NDUFB7_HF_{2026_04_20}"
DATA_RAW="$PROJECT_ROOT/01_data/01_raw_geo"
LOG_DIR="$PROJECT_ROOT/05_logs"
mkdir -p "$DATA_RAW" "$LOG_DIR"

MASTER_LOG="$LOG_DIR/batch_download_$(date +%Y%m%d_%H%M%S).log"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] === 批量GEO下载启动 ===" | tee -a "$MASTER_LOG"

# GEO FTP路径计算（修复版：取GSE数字部分前3位）
get_ftp_dir() {
    local num=${1#GSE}
    local prefix=${num:0:3}
    echo "GSE${prefix}nnn"
}

# 下载函数
dl_gse() {
    local gse=$1
    local priority=$2
    local target="$DATA_RAW/$gse"
    mkdir -p "$target"
    
    local ftp_dir=$(get_ftp_dir "$gse")
    local ftp_url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/$ftp_dir/$gse"
    
    echo "" | tee -a "$MASTER_LOG"
    echo "[$(date '+%H:%M:%S')] [$priority] 开始: $gse" | tee -a "$MASTER_LOG"
    
    # 1. 下载Series Matrix（断点续传）
    local sm_file="$target/${gse}_series_matrix.txt.gz"
    if [ -f "$sm_file" ] && gunzip -t "$sm_file" 2>/dev/null; then
        echo "  -> Series Matrix已存在且完好，跳过" | tee -a "$MASTER_LOG"
    else
        echo "  -> 下载Series Matrix..." | tee -a "$MASTER_LOG"
        wget -c -t 0 --timeout=120 --progress=dot:giga \
            -O "$sm_file.tmp" \
            "$ftp_url/matrix/${gse}_series_matrix.txt.gz" \
            2>>"$LOG_DIR/${gse}_wget.log" && \
        mv "$sm_file.tmp" "$sm_file" && \
        echo "  -> Series Matrix完成" | tee -a "$MASTER_LOG" || \
        echo "  -> Series Matrix失败（可能无matrix或网络问题）" | tee -a "$MASTER_LOG"
    fi
    
    # 2. 检查并下载Supplementary（仅对单细胞/大数据集）
    local suppl_url="$ftp_url/suppl/"
    if wget --spider --timeout=60 "$suppl_url" 2>/dev/null; then
        echo "  -> 发现Supplementary目录，下载关键文件..." | tee -a "$MASTER_LOG"
        mkdir -p "$target/suppl"
        # 只下载.gz文件，避免递归太深
        wget -c -t 0 --timeout=120 --progress=dot:giga -r -np -nH --cut-dirs=3 -A "*.gz" \
            -P "$target/suppl" \
            "$suppl_url" \
            2>>"$LOG_DIR/${gse}_suppl.log" || true
    fi
    
    # 3. 记录状态
    local fcount=$(ls -1 "$target" 2>/dev/null | wc -l)
    local dsize=$(du -sh "$target" 2>/dev/null | cut -f1)
    echo "  -> 完成: $dsize ($fcount files)" | tee -a "$MASTER_LOG"
    echo "$gse,$priority,$(date '+%H:%M:%S'),$dsize,$fcount" >> "$LOG_DIR/download_status.csv"
    
    sleep 3
}

# GSE列表（与你的研究方案对齐）
P0="GSE141910 GSE57338 GSE168742"
P1="GSE183852 GSE116250 GSE46224 GSE5406 GSE271946"
P2="GSE315590 GSE217494 GSE270788"

echo "GSE,Priority,Time,Size,Files" > "$LOG_DIR/download_status.csv"

# 执行下载
for gse in $P0; do dl_gse "$gse" "P0_must"; done
for gse in $P1; do dl_gse "$gse" "P1_urgent"; done
for gse in $P2; do dl_gse "$gse" "P2_extended"; done

echo "" | tee -a "$MASTER_LOG"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] === 批量下载任务结束 ===" | tee -a "$MASTER_LOG"
echo "[日志] $MASTER_LOG"
echo "[状态] $LOG_DIR/download_status.csv"
