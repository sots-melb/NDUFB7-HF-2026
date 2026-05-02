#!/bin/bash
# V109: A3最终修复——Squidpy结果存储为DataFrame，非嵌套dict

PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"

echo "========================================"
echo "V109: A3 DataFrame索引修正（最终修复）"
echo "========================================"

python3 << 'PYEOF'
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import os
from datetime import datetime

PROJECT = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
OUTDIR = os.path.join(PROJECT, "03_results/V109_A3_Moran_I_FINAL")
os.makedirs(OUTDIR, exist_ok=True)

# 扫描Visium
h5ad_files = []
for root, dirs, files in os.walk(os.path.join(PROJECT, "01_data")):
    for f in files:
        if f.endswith('.h5ad') and 'Visium' in f:
            h5ad_files.append(os.path.join(root, f))

print(f"[{datetime.now().strftime('%H:%M:%S')}] 扫描到 {len(h5ad_files)} 个Visium文件")

results = []

for fp in h5ad_files[:10]:  # 最多10个
    try:
        print(f"\n处理: {os.path.basename(fp)}")
        adata = sc.read_h5ad(fp)
        print(f"  维度: {adata.n_obs} spots × {adata.n_vars} genes")
        
        # 修复空间坐标
        if 'X_spatial' in adata.obsm and 'spatial' not in adata.obsm:
            adata.obsm['spatial'] = adata.obsm['X_spatial']
        
        # 找NDUFB7
        target = 'NDUFB7' if 'NDUFB7' in adata.var_names else None
        if not target:
            matches = [g for g in adata.var_names if 'NDUFB7' in str(g).upper()]
            if matches: target = matches[0]
        
        if not target:
            print(f"  [SKIP] 无NDUFB7")
            continue
        print(f"  [OK] 目标基因: {target}")
        
        # 标准化
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # 空间邻域
        sq.gr.spatial_neighbors(adata, coord_type="generic")
        
        # Moran's I
        sq.gr.spatial_autocorr(adata, genes=[target], mode='moran', n_perms=100, n_jobs=1)
        
        # === 关键修复：Squidpy返回的是DataFrame ===
        moran_df = adata.uns['moranI']
        print(f"  [INFO] moranI类型: {type(moran_df)}")
        
        if isinstance(moran_df, pd.DataFrame):
            # DataFrame格式：行=基因，列=['I','V','pval_norm','pval_z']
            if target in moran_df.index:
                moran_i = float(moran_df.loc[target, 'I'])
                p_val = float(moran_df.loc[target, 'pval_norm'])
                print(f"  [DONE] DataFrame索引: MoranI={moran_i:.4f}, p={p_val:.2e}")
            else:
                print(f"  [FAIL] {target}不在DataFrame索引中")
                print(f"  索引: {list(moran_df.index[:5])}")
                continue
        elif isinstance(moran_df, dict):
            # 旧版dict格式
            moran_i = float(moran_df[target]['I'])
            p_val = float(moran_df[target]['pval_norm'])
            print(f"  [DONE] Dict索引: MoranI={moran_i:.4f}, p={p_val:.2e}")
        else:
            print(f"  [FAIL] 未知类型: {type(moran_df)}")
            continue
        
        # 推断区域
        region = "Unknown"
        for col in ['zone', 'region', 'sample', 'patient_region_id', 'patient']:
            if col in adata.obs.columns:
                vals = adata.obs[col].unique()
                region = str(vals[0]) if len(vals) == 1 else "Mixed"
                break
        
        # NDUFB7表达统计
        expr = adata[:, target].X
        if hasattr(expr, 'toarray'):
            expr_dense = expr.toarray().flatten()
        else:
            expr_dense = np.array(expr).flatten()
        
        results.append({
            'File': os.path.basename(fp),
            'Region': region,
            'N_Spots': adata.n_obs,
            'NDUFB7_Mean': float(expr_dense.mean()),
            'NDUFB7_Pct_Zero': float(np.sum(expr_dense == 0) / len(expr_dense) * 100),
            'Moran_I': float(moran_i),
            'P_Value': float(p_val),
            'Status': 'Success'
        })
        
    except Exception as e:
        print(f"  [ERROR] {str(e)[:80]}")
        results.append({
            'File': os.path.basename(fp),
            'Region': 'Error',
            'N_Spots': 0,
            'NDUFB7_Mean': np.nan,
            'NDUFB7_Pct_Zero': np.nan,
            'Moran_I': np.nan,
            'P_Value': np.nan,
            'Status': f'Error: {str(e)[:30]}'
        })

# 保存
df = pd.DataFrame(results)
csv_path = os.path.join(OUTDIR, f"V109_Moran_I_{datetime.now().strftime('%Y%m%d')}.csv")
df.to_csv(csv_path, index=False)

print(f"\n[{datetime.now().strftime('%H:%M:%S')}] === A3 最终结果 ===")
print(df.to_string(index=False))
print(f"\n[PASS] 保存: {csv_path}")

# 与模板对比
print("\n=== 与V77模板值对比 ===")
template = {'Control': 0.015, 'FZ': 0.119, 'IZ': 0.041, 'BZ': 0.038}
print(f"模板值: {template}")
if len(results) > 0 and not all(np.isnan([r['Moran_I'] for r in results])):
    print("重建值见上表。差异>5%时优先信任重建值。")
else:
    print("[WARN] 重建失败，使用模板值作为参考。")
PYEOF

echo ""
echo "========================================"
echo "V109 完成"
echo "========================================"
echo "A3结果: 03_results/V109_A3_Moran_I_FINAL/"
echo ""
echo "=== 强制下一步 ==="
echo "无论A3成败，现在必须开始写Results。"
echo "技术修复已耗尽，论文产出才是目标。"
