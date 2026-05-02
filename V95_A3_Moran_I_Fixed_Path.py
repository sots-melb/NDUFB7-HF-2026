#!/usr/bin/env python3
"""
V95_A3: Moran's I 重建（固定路径）
从Project内Visium h5ad计算空间自相关
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import os
from datetime import datetime

PROJECT_DIR = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "03_results/A3_Moran_I")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 严格固定路径
VISIUM_FILES = [
    "01_data/02_spatial/Visium_RZ_BZ_P12.h5ad",
    "01_data/02_spatial/Visium_IZ_BZ_P2.h5ad",
    "01_data/02_spatial/Visium_Control_P1.h5ad"
]

print(f"[{datetime.now().strftime('%H:%M:%S')}] === V95 A3: Moran's I 重建 ===")
print(f"输出目录: {OUTPUT_DIR}")

results = []

for f in VISIUM_FILES:
    path = os.path.join(PROJECT_DIR, f)
    if not os.path.exists(path):
        print(f"[SKIP] 不存在: {f}")
        continue
    
    print(f"\n[{datetime.now().strftime('%H:%M:%S')}] 处理: {os.path.basename(path)}")
    
    try:
        adata = sc.read_h5ad(path)
        print(f"  维度: {adata.n_obs} spots × {adata.n_vars} genes")
        
        # 质控
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        
        # 标准化
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # 检查NDUFB7
        if "NDUFB7" not in adata.var_names:
            print(f"  [WARN] NDUFB7 不在基因列表中")
            continue
        
        # 空间邻域图
        sq.gr.spatial_neighbors(adata)
        
        # Moran's I
        sq.gr.spatial_autocorr(adata, genes=["NDUFB7"], mode="moran")
        
        moran_i = adata.uns["moranI"]["NDUFB7"]["I"]
        p_val = adata.uns["moranI"]["NDUFB7"]["pval_norm"]
        
        # 推断区域
        region = "Unknown"
        if "zone" in adata.obs:
            region = str(adata.obs["zone"].iloc[0])
        elif "region" in adata.obs:
            region = str(adata.obs["region"].iloc[0])
        elif "sample" in adata.obs:
            region = str(adata.obs["sample"].iloc[0])
        
        results.append({
            "File": os.path.basename(path),
            "Region": region,
            "N_Spots": adata.n_obs,
            "NDUFB7_Mean": float(adata[:, "NDUFB7"].X.mean()),
            "NDUFB7_Pct": float((adata[:, "NDUFB7"].X > 0).mean() * 100),
            "Moran_I": float(moran_i),
            "P_Value": float(p_val),
            "Status": "Success"
        })
        
        print(f"  [DONE] Moran's I = {moran_i:.4f}, p = {p_val:.2e}")
        
    except Exception as e:
        print(f"  [ERROR] {str(e)}")
        results.append({
            "File": os.path.basename(path),
            "Region": "Error",
            "N_Spots": 0,
            "NDUFB7_Mean": np.nan,
            "NDUFB7_Pct": np.nan,
            "Moran_I": np.nan,
            "P_Value": np.nan,
            "Status": f"Error: {str(e)[:50]}"
        })

# 保存结果
df = pd.DataFrame(results)
csv_path = os.path.join(OUTPUT_DIR, f"V95_Moran_I_{datetime.now().strftime('%Y%m%d')}.csv")
df.to_csv(csv_path, index=False)

print(f"\n[{datetime.now().strftime('%H:%M:%S')}] === 结果汇总 ===")
print(df.to_string(index=False))
print(f"\n[PASS] 结果保存: {csv_path}")

# 与模板值对比
print("\n=== 与V77模板值对比 ===")
template = {"Control": 0.015, "FZ": 0.119, "IZ": 0.041, "BZ": 0.038}
print("模板值:", template)
print("重建值见上表。差异>5%时优先信任重建值。")
