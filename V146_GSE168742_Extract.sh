#!/bin/bash
PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"
OUTDIR="$PROJECT/03_results/V146_GSE168742_NDUFB7"
mkdir -p "$OUTDIR"
LOG="$OUTDIR/extraction.log"

echo "=== V146: GSE168742 NDUFB7提取 ===" | tee "$LOG"
echo "Start: $(date)" | tee -a "$LOG"

RAW_TAR="$HOME/Downloads/GSE168742_RAW.tar"
if [ ! -f "$RAW_TAR" ]; then
  echo "[ERROR] GSE168742_RAW.tar not found in Downloads" | tee -a "$LOG"
  exit 1
fi

# 解压到临时目录
TMPDIR="$OUTDIR/tmp_extract"
mkdir -p "$TMPDIR"
echo "[1/4] Extracting RAW.tar..." | tee -a "$LOG"
tar -xf "$RAW_TAR" -C "$TMPDIR" 2>&1 | head -20 | tee -a "$LOG"

# 检查内容
echo "[2/4] Checking extracted files..." | tee -a "$LOG"
ls -lh "$TMPDIR" | tee -a "$LOG"

# 尝试用Python读取并提取NDUFB7
echo "[3/4] Extracting NDUFB7..." | tee -a "$LOG"
python3 << 'PYEOF' 2>&1 | tee -a "$LOG"
import os, sys, glob, warnings
warnings.filterwarnings("ignore")
import pandas as pd

tmp = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/03_results/V146_GSE168742_NDUFB7/tmp_extract")
out = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/03_results/V146_GSE168742_NDUFB7")

# 寻找可能的表达矩阵文件
candidates = []
for ext in ["*.csv", "*.txt", "*.tsv", "*.mtx", "*.h5", "*.h5ad"]:
  candidates.extend(glob.glob(os.path.join(tmp, "**", ext), recursive=True))

print(f"Found {len(candidates)} candidate files")
for c in candidates[:10]:
  print("  ", os.path.basename(c), os.path.getsize(c))

# 简化策略：如果有h5ad直接读取；如果有10x格式读取；如果有txt/csv读取
ndufb7_vals = []
sample_ids = []

if candidates:
    # 尝试读取第一个csv/txt
    for c in candidates:
        if c.endswith(('.csv','.txt','.tsv')) and os.path.getsize(c) > 1000:
            try:
                sep = '\t' if c.endswith('.tsv') else (',' if c.endswith('.csv') else '\t')
                df = pd.read_csv(c, sep=sep, index_col=0, nrows=5)
                if 'NDUFB7' in df.index or 'NDUFB7' in str(df.index):
                    df_full = pd.read_csv(c, sep=sep, index_col=0)
                    if 'NDUFB7' in df_full.index:
                        ndufb7_vals = df_full.loc['NDUFB7'].values.flatten().tolist()
                        sample_ids = df_full.columns.tolist()
                        print(f"[PASS] NDUFB7 extracted from {os.path.basename(c)}: {len(ndufb7_vals)} samples")
                        break
                    # 尝试大小写
                    idx_match = [i for i in df_full.index if 'NDUFB7' in str(i).upper()]
                    if idx_match:
                        ndufb7_vals = df_full.loc[idx_match[0]].values.flatten().tolist()
                        sample_ids = df_full.columns.tolist()
                        print(f"[PASS] NDUFB7 extracted (fuzzy match {idx_match[0]}): {len(ndufb7_vals)} samples")
                        break
            except Exception as e:
                print(f"  [SKIP] {os.path.basename(c)}: {str(e)[:50]}")
                continue

if ndufb7_vals:
    out_df = pd.DataFrame({'Sample': sample_ids, 'NDUFB7': ndufb7_vals})
    out_df.to_csv(os.path.join(out, "V146_GSE168742_NDUFB7.csv"), index=False)
    print(f"[DONE] Saved V146_GSE168742_NDUFB7.csv, n={len(ndufb7_vals)}")
else:
    print("[FAIL] Could not extract NDUFB7 automatically. Manual inspection needed.")
    # 保存文件列表供手动检查
    with open(os.path.join(out, "file_manifest.txt"), 'w') as f:
        for c in candidates:
            f.write(f"{os.path.getsize(c)}\t{c}\n")
    print("[INFO] File manifest saved for manual inspection")
PYEOF

echo "[4/4] Cleanup..." | tee -a "$LOG"
rm -rf "$TMPDIR"
echo "End: $(date)" | tee -a "$LOG"
