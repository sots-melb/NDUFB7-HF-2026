#!/bin/bash
PROJECT_ROOT="$HOME/Projects/NDUFB7_HF_{2026_04_20}"
DATA_RAW="$PROJECT_ROOT/01_data/01_raw_geo"
LOG_DIR="$PROJECT_ROOT/05_logs"
mkdir -p "$DATA_RAW" "$LOG_DIR"

LOG="$LOG_DIR/suppl_download_$(date +%Y%m%d_%H%M%S).log"
echo "[$(date '+%H:%M:%S')] 开始补下载Supplementary..." | tee -a "$LOG"

# 只下载有processed counts的GSE（从GEO页面确认过有suppl文件）
for gse in GSE116250 GSE46224 GSE315590; do
    target="$DATA_RAW/$gse"
    mkdir -p "$target"
    
    # 尝试下载常见的processed文件名
    for pattern in "${gse}_gene_counts.csv.gz" "${gse}_counts.txt.gz" "${gse}_fpkm.txt.gz" "${gse}_tpm.txt.gz" "${gse}_RAW.tar"; do
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/${gse:0:6}nnn/$gse/suppl/$pattern"
        echo "[$gse] 尝试: $pattern" | tee -a "$LOG"
        wget -c -t 2 --timeout=60 --progress=dot:giga -O "$target/$pattern" "$url" 2>/dev/null && \
        echo "[$gse] ✅ $pattern 成功" | tee -a "$LOG" || \
        echo "[$gse] ❌ $pattern 失败" | tee -a "$LOG"
    done
    sleep 2
done

echo "[$(date '+%H:%M:%S')] 结束" | tee -a "$LOG"
