import scanpy as sc
import pandas as pd
import os, glob, warnings
warnings.filterwarnings("ignore")

PROJECT = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
spatial_dir = os.path.join(PROJECT, "01_data/02_spatial")
outdir = os.path.join(PROJECT, "03_results/V147_Fig2_Publication")
os.makedirs(outdir, exist_ok=True)

# 找3个Visium文件（control优先）
h5ad_files = glob.glob(os.path.join(spatial_dir, "**/*.h5ad"), recursive=True)
h5ad_files = [f for f in h5ad_files if 'control' in f.lower() or 'ctrl' in f.lower()][:3]
if len(h5ad_files) < 3:
    h5ad_files = glob.glob(os.path.join(spatial_dir, "**/*.h5ad"), recursive=True)[:3]

print(f"[INFO] Processing {len(h5ad_files)} Visium files")

all_data = []
for f in h5ad_files:
    try:
        ad = sc.read_h5ad(f)
        has_ndufb7 = "NDUFB7" in ad.var_names
        has_umap = "X_umap" in ad.obsm
        
        if has_umap:
            coords = ad.obsm["X_umap"][:, :2]
        elif "X_pca" in ad.obsm:
            coords = ad.obsm["X_pca"][:, :2]
        else:
            continue
            
        y = None
        if has_ndufb7:
            y = ad[:, "NDUFB7"].X.toarray().flatten() if hasattr(ad[:, "NDUFB7"].X, "toarray") else ad[:, "NDUFB7"].X.flatten()
        
        df = pd.DataFrame({
            'File': os.path.basename(f),
            'UMAP1': coords[:,0],
            'UMAP2': coords[:,1],
            'NDUFB7': y if y is not None else 0,
            'Has_NDUFB7': has_ndufb7
        })
        all_data.append(df)
        print(f"  [PASS] {os.path.basename(f)}: n={ad.n_obs}, NDUFB7={has_ndufb7}")
    except Exception as e:
        print(f"  [SKIP] {os.path.basename(f)}: {str(e)[:50]}")

if all_data:
    combined = pd.concat(all_data, ignore_index=True)
    combined.to_csv(os.path.join(outdir, "V147A_visium_umap_ndufb7.csv"), index=False)
    print(f"[DONE] Saved V147A_visium_umap_ndufb7.csv, total spots={len(combined)}")
else:
    print("[FAIL] No Visium data extracted")
