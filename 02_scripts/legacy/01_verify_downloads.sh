#!/bin/bash
PROJECT_ROOT="$HOME/Projects/NDUFB7_HF_{2026_04_20}"
DATA_RAW="$PROJECT_ROOT/01_data/01_raw_geo"
REPORT="$PROJECT_ROOT/03_results/01_qc_reports/download_verification_$(date +%Y%m%d_%H%M%S).txt"

echo "=== GEO数据完整性验证报告 ===" | tee "$REPORT"
echo "生成时间: $(date)" | tee -a "$REPORT"
echo "" | tee -a "$REPORT"

for d in GSE141910 GSE57338 GSE168742 GSE183852 GSE116250 GSE46224 GSE5406 GSE271946 GSE315590 GSE217494 GSE270788; do
    dir_path="$DATA_RAW/$d"
    if [ ! -d "$dir_path" ]; then
        echo "[$d] ❌ 目录不存在" | tee -a "$REPORT"
        continue
    fi
    
    # 查找最大文件
    largest=$(find "$dir_path" -type f -printf '%s %p\n' 2>/dev/null | sort -nr | head -1)
    if [ -z "$largest" ]; then
        echo "[$d] ⚠️  目录为空" | tee -a "$REPORT"
        continue
    fi
    
    size_b=$(echo "$largest" | awk '{print $1}')
    size_mb=$(awk "BEGIN {printf \"%.1f\", $size_b/1024/1024}")
    fname=$(echo "$largest" | awk '{$1=""; print $0}' | sed 's/^ //' | xargs basename)
    
    # gzip校验
    if gunzip -t "$(echo "$largest" | awk '{$1=""; print $0}' | sed 's/^ //')" 2>/dev/null; then
        gz_status="gzip校验通过"
    else
        gz_status="gzip校验失败"
    fi
    
    # 大小判断
    size_int=$(echo "$size_mb" | cut -d. -f1)
    if [ "$size_int" -gt 5 ] 2>/dev/null; then
        status="✅ 正常"
    elif [ "$size_int" -gt 0 ] 2>/dev/null; then
        status="⚠️  偏小（可能仅metadata）"
    else
        status="⚠️  极小"
    fi
    
    echo "[$d] $status | ${size_mb}MB | $fname | $gz_status" | tee -a "$REPORT"
done

echo "" | tee -a "$REPORT"
echo "[完成] 报告保存: $REPORT" | tee -a "$REPORT"
