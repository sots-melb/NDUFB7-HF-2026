#!/usr/bin/env python3
import scanpy as sc
import pandas as pd
import os, tarfile, shutil, warnings
warnings.filterwarnings("ignore")

PROJECT = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
outdir = os.path.join(PROJECT, "03_results/V148_GSE168742_Fixed")
os.makedirs(outdir, exist_ok=True)

RAW_TAR = os.path.expanduser("~/Downloads/GSE168742_RAW.tar")
LOG = os.path.join(outdir, "extraction.log")

def log(msg):
    print(msg)
    with open(LOG, 'a') as f:
        f.write(msg + "\n")

log("=== V148: GSE168742 10x格式NDUFB7提取 ===")

if not os.path.exists(RAW_TAR):
    log("[ERROR] GSE168742_RAW.tar not found")
    exit(1)

tmpdir = os.path.join(outdir, "tmp_extract")
if os.path.exists(tmpdir):
    shutil.rmtree(tmpdir)
os.makedirs(tmpdir, exist_ok=True)

log("[1/3] Extracting tar...")
with tarfile.open(RAW_TAR, 'r') as tar:
    tar.extractall(path=tmpdir)
log("  Extracted")

# 关键修复：递归寻找10x目录（嵌套在子目录中）
log("[2/3] Scanning for 10x directories...")
tenx_dirs = []
for root, dirs, files in os.walk(tmpdir):
    has_barcode = any(f.endswith('barcodes.tsv.gz') for f in files)
    has_features = any(f.endswith('features.tsv.gz') for f in files) or any(f.endswith('genes.tsv.gz') for f in files)
    has_matrix = any(f.endswith('matrix.mtx.gz') for f in files)
    if has_barcode and (has_features or has_matrix):
        tenx_dirs.append(root)
        log(f"  Found 10x dir: {root} (files: {len(files)})")

if not tenx_dirs:
    log("[FAIL] No 10x directories found")
    all_files = []
    for root, dirs, files in os.walk(tmpdir):
        for f in files:
            all_files.append(os.path.join(root, f))
    with open(os.path.join(outdir, "file_manifest.txt"), 'w') as f:
        for af in sorted(all_files):
            f.write(f"{os.path.getsize(af)}\t{af}\n")
    log("[INFO] File manifest saved")
    shutil.rmtree(tmpdir)
    exit(1)

# 读取所有10x样本
log("[3/3] Reading 10x data and extracting NDUFB7...")
all_ndufb7 = []

for td in tenx_dirs[:5]:  # 最多5个样本
    try:
        log(f"  Reading {os.path.basename(td)}...")
        ad = sc.read_10x_mtx(td, var_names='gene_symbols', cache=False)
        log(f"    Shape: {ad.n_obs} cells × {ad.n_vars} genes")
        
        if 'NDUFB7' in ad.var_names:
            y = ad[:, 'NDUFB7'].X.toarray().flatten() if hasattr(ad[:, 'NDUFB7'].X, 'toarray') else ad[:, 'NDUFB7'].X.flatten()
            all_ndufb7.append({
                'Sample': os.path.basename(td),
                'N_Cells': ad.n_obs,
                'NDUFB7_Mean': float(y.mean()),
                'NDUFB7_Median': float(pd.Series(y).median()),
                'NDUFB7_ZeroPct': float((y == 0).mean() * 100)
            })
            log(f"    [PASS] NDUFB7 mean={y.mean():.3f}, zero%={(y==0).mean()*100:.1f}%")
        else:
            # 尝试ensembl ID
            ensembl_match = [v for v in ad.var_names if 'NDUFB7' in str(v).upper()]
            if ensembl_match:
                y = ad[:, ensembl_match[0]].X.toarray().flatten() if hasattr(ad[:, ensembl_match[0]].X, 'toarray') else ad[:, ensembl_match[0]].X.flatten()
                log(f"    [PASS] NDUFB7 matched via {ensembl_match[0]}")
            else:
                log(f"    [WARN] NDUFB7 not found")
                all_ndufb7.append({
                    'Sample': os.path.basename(td), 'N_Cells': ad.n_obs,
                    'NDUFB7_Mean': 0, 'NDUFB7_Median': 0, 'NDUFB7_ZeroPct': 100
                })
    except Exception as e:
        log(f"    [ERROR] {str(e)[:100]}")

if all_ndufb7:
    df = pd.DataFrame(all_ndufb7)
    df.to_csv(os.path.join(outdir, "V148_GSE168742_NDUFB7_summary.csv"), index=False)
    log(f"[PASS] Saved summary for {len(df)} samples")
else:
    log("[FAIL] No NDUFB7 extracted")

shutil.rmtree(tmpdir)
log(f"[DONE] Output: {outdir}")
