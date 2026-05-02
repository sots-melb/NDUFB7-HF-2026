#!/bin/bash
OUTDIR="$PROJECT/03_results/V168_PanDeath_Solidity/V168A_DOX"
mkdir -p "$OUTDIR"
echo "========================================"
echo "V168A: GSE157282 DOX心肌毒性 — NDUFB7诊断"
echo "========================================"

TAR="$HOME/Downloads/GSE157282_RAW.tar"
SM="$HOME/Downloads/GSE157282_series_matrix.txt.gz"

# 1. 解析series matrix分组
echo ""
echo "[1/3] Series matrix分组诊断..."
if [ -f "$SM" ]; then
  zcat "$SM" | grep -E "^!Sample_title|^!Sample_source_name_ch1|^!Sample_characteristics" | head -20
  echo "---"
  # 提取DOX vs Control标签
  zcat "$SM" | grep -i "dox\\|control\\|vehicle\\|saline" | head -10
else
  echo "[FAIL] Series matrix not found"
fi

# 2. 检查tar内容
echo ""
echo "[2/3] RAW tar内容诊断..."
if [ -f "$TAR" ]; then
  tar -tf "$TAR" | head -20
  echo "---"
  echo "文件总数: $(tar -tf "$TAR" | wc -l)"
  
  # 寻找表达矩阵或count文件
  echo ""
  echo "寻找表达/count文件..."
  tar -tf "$TAR" | grep -iE "count|expr|matrix|txt|csv" | head -10
else
  echo "[FAIL] RAW tar not found"
fi

# 3. 尝试直接提取NDUFB7（如果有featureCounts或表达矩阵）
echo ""
echo "[3/3] 尝试提取NDUFB7..."
# GSE157282通常是小鼠RNA-seq，可能包含featureCounts
TMP=$(mktemp -d)
tar -xf "$TAR" -C "$TMP" 2>/dev/null
FC=$(find "$TMP" -name "*featureCounts*" -o -name "*counts.txt*" | head -1)
if [ -n "$FC" ]; then
  echo "Found count file: $FC"
  head -5 "$FC"
  echo "---"
  grep -i "ndufb7\|Ndufb7" "$FC" | head -3 || echo "NDUFB7 not found in counts"
else
  echo "[INFO] No featureCounts found — may need alignment from FASTQ"
fi
rm -rf "$TMP"

echo ""
echo "[DONE] V168A诊断完成"
echo "  下一步: 如果有表达矩阵，运行R提取脚本"
