#!/bin/bash
OUTDIR="$HOME/Projects/NDUFB7_HF_2026_04_20/03_results/V148_GSE168742_Fixed"
TMPDIR="$OUTDIR/tmp_extract"
mkdir -p "$TMPDIR"

echo "=== V148诊断：查看tar内容和目录结构 ==="
cd ~/Downloads && tar -tf GSE168742_RAW.tar | head -30
echo "---"
tar -xf GSE168742_RAW.tar -C "$TMPDIR"
find "$TMPDIR" -type f | head -20
echo "---"
find "$TMPDIR" -name "*.mtx.gz" -o -name "*.tsv.gz" | head -10
