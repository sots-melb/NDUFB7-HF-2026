#!/bin/bash
PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"
GEO="$PROJECT/01_data/01_raw_geo"
LOG="$PROJECT/04_logs"
mkdir -p "$LOG"
REPORT="$LOG/data_verification_report.txt"
echo "═══════════════════════════════════════════════════════════════" > "$REPORT"
echo "  NDUFB7_Mito_2026 - 数据校验报告" >> "$REPORT"
echo "  时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$REPORT"
echo "═══════════════════════════════════════════════════════════════" >> "$REPORT"
echo "" >> "$REPORT"
TOTAL=0; OK=0; FAIL=0
for gse_dir in $(ls -1 "$GEO" 2>/dev/null | sort); do
    TOTAL=$((TOTAL + 1))
    local matrix_file="$GEO/$gse_dir/${gse_dir}_series_matrix.txt.gz"
    if [ -f "$matrix_file" ] && [ -s "$matrix_file" ] && gunzip -t "$matrix_file" 2>/dev/null; then
        local sz=$(du -sh "$matrix_file" | cut -f1)
        echo "  ✅ $gse_dir ($sz)" | tee -a "$REPORT"
        OK=$((OK + 1))
    else
        echo "  ❌ $gse_dir" | tee -a "$REPORT"
        FAIL=$((FAIL + 1))
    fi
done
echo "" | tee -a "$REPORT"
echo "总计: $TOTAL | ✅成功:$OK | ❌失败:$FAIL" | tee -a "$REPORT"
echo "磁盘: $(df -h $PROJECT | awk 'NR==2 {print $4}') 可用" | tee -a "$REPORT"
echo "报告: $REPORT"
