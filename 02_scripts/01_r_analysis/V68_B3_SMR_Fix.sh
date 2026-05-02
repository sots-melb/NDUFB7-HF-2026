#!/bin/bash
SMR="$HOME/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/smr_software/smr"
OUT="$HOME/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables"
DATA="$HOME/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data"

echo "=== 检查数据文件 ==="
ls -lh "$DATA/smr_input/HERMES_HF_SMR.ma.gz" 2>/dev/null || echo "HERMES .ma.gz 未找到"
ls -lh "$DATA/hermes_v1/HERMES_v1_for_SMR.ma.gz" 2>/dev/null || echo "HERMES_v1 .ma.gz 未找到"

# 解压HERMES
if [ -f "$DATA/smr_input/HERMES_HF_SMR.ma.gz" ]; then
  echo ">>> 解压 HERMES..."
  gunzip -c "$DATA/smr_input/HERMES_HF_SMR.ma.gz" > "$DATA/smr_input/HERMES_HF_SMR.ma" 2>/dev/null && echo "✅ 解压成功" || echo "❌ 解压失败"
elif [ -f "$DATA/hermes_v1/HERMES_v1_for_SMR.ma.gz" ]; then
  echo ">>> 解压 HERMES_v1..."
  gunzip -c "$DATA/hermes_v1/HERMES_v1_for_SMR.ma.gz" > "$DATA/smr_input/HERMES_HF_SMR.ma" 2>/dev/null && echo "✅ 解压成功" || echo "❌ 解压失败"
fi

# 检查eQTLGen格式
echo ">>> eQTLGen提取文件预览:"
head -3 "$DATA/eqtlgen_NDUFB7_extracted.csv" 2>/dev/null || echo "  不存在"

# 检查GTEx格式
echo ">>> GTEx预览:"
head -3 "$DATA/gtex_NDUFB7_HLV_v11.csv" 2>/dev/null || echo "  不存在"

echo "=== SMR数据审计完成 ==="
