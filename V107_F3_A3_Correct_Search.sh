#!/bin/bash
# V107: F3用gene_name列 + A3全列表搜索NDUFB7

PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"

echo "========================================"
echo "V107: 最终修复（正确列/搜索范围）"
echo "========================================"

# ========================================
# F3: 用gene_name列提取NDUFB7
# ========================================
echo ""
echo ">>> [F3] 从featureCounts的gene_name列提取"

python3 << 'PYEOF'
import pandas as pd
import os
from scipy import stats

EXTRACT_DIR = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE243655")
OUTDIR = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/03_results/V107_F3_Fer1")
os.makedirs(OUTDIR, exist_ok=True)

groups = {
    'GSM7792365': 'DMSO', 'GSM7792366': 'DMSO', 
    'GSM7792367': 'DMSO', 'GSM7792368': 'DMSO',
    'GSM7792369': 'Fer-1', 'GSM7792370': 'Fer-1',
    'GSM7792371': 'Fer-1', 'GSM7792372': 'Fer-1'
}

results = []
for gsm, grp in groups.items():
    import glob
    pattern = os.path.join(EXTRACT_DIR, f"{gsm}*featureCounts.txt.gz")
    files = glob.glob(pattern)
    if not files:
        continue
    
    # featureCounts格式：Geneid Chr Start End Strand Length gene_name sample_counts
    df = pd.read_csv(files[0], sep='\t', skiprows=1)
    
    # 找gene_name列（倒数第二列）
    cols = df.columns.tolist()
    gene_name_col = [c for c in cols if 'gene_name' in c.lower()][0]
    counts_col = cols[-1]  # 最后一列是counts
    
    # 在gene_name列搜索NDUFB7
    match = df[df[gene_name_col].str.contains('^NDUFB7$', case=False, na=False, regex=True)]
    
    if match.empty:
        # 尝试模糊匹配
        match = df[df[gene_name_col].str.contains('NDUFB7', case=False, na=False)]
        if not match.empty:
            print(f"  [WARN] {gsm}: 模糊匹配到 {match.iloc[0][gene_name_col]}")
    
    if not match.empty:
        counts = int(match.iloc[0][counts_col])
        gene = match.iloc[0][gene_name_col]
        results.append({'GSM': gsm, 'Group': grp, 'Gene': gene, 'Counts': counts})
        print(f"  [OK] {gsm} ({grp}): {gene} = {counts}")
    else:
        # 诊断：列出所有NDUF开头的基因
        nduf_genes = df[df[gene_name_col].str.contains('^NDUF', case=False, na=False, regex=True)][gene_name_col].unique()
        print(f"  [FAIL] {gsm}: NDUFB7未找到。可用NDUF基因: {', '.join(nduf_genes[:10])}")

print(f"\n提取完成: {len(results)}/8 样本")

if len(results) == 8:
    df_res = pd.DataFrame(results)
    dmso = df_res[df_res['Group']=='DMSO']['Counts'].values
    fer1 = df_res[df_res['Group']=='Fer-1']['Counts'].values
    
    t_stat, p_val = stats.ttest_rel(dmso, fer1)
    delta = fer1.mean() - dmso.mean()
    pct = delta / dmso.mean() * 100
    
    print(f"\n=== F3 配对t-test ===")
    print(f"DMSO: {dmso.mean():.1f} ± {dmso.std():.1f}")
    print(f"Fer-1: {fer1.mean():.1f} ± {fer1.std():.1f}")
    print(f"Δ = {delta:.1f} ({pct:+.1f}%)")
    print(f"t(3) = {t_stat:.3f}")
    print(f"p = {p_val:.4f}")
    
    verdict = "PASS" if (p_val < 0.05 and fer1.mean() > dmso.mean()) else ("PARTIAL" if p_val < 0.05 else "FAIL")
    print(f"\n[{verdict}] Fer-1 {'上调' if fer1.mean() > dmso.mean() else '下调'}NDUFB7")
    
    df_res.to_csv(os.path.join(OUTDIR, "V107_Fer1_NDUFB7.csv"), index=False)
    pd.DataFrame({
        'Patient': [1,2,3,4], 'DMSO': dmso, 'Fer1': fer1, 'Delta': fer1-dmso
    }).to_csv(os.path.join(OUTDIR, "V107_Fer1_paired.csv"), index=False)
    print(f"[DONE] 保存: {OUTDIR}")
else:
    print(f"[FAIL] 样本不足，无法统计")
PYEOF

# ========================================
# A3: 在h5ad全列表中搜索NDUFB7
# ========================================
echo ""
echo ">>> [A3] h5ad全列表搜索NDUFB7"

python3 << 'PYEOF'
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import os

PROJECT = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
OUTDIR = os.path.join(PROJECT, "03_results/V107_A3_Moran_I")
os.makedirs(OUTDIR, exist_ok=True)

# 扫描Visium文件
h5ad_files = []
for root, dirs, files in os.walk(os.path.join(PROJECT, "01_data")):
    for f in files:
        if f.endswith('.h5ad') and 'Visium' in f:
            h5ad_files.append(os.path.join(root, f))

print(f"扫描到 {len(h5ad_files)} 个Visium文件")

# 先检查第一个文件是否有NDUFB7
if h5ad_files:
    adata = sc.read_h5ad(h5ad_files[0])
    all_genes = adata.var_names.tolist()
    
    # 全列表搜索
    ndufb7_matches = [g for g in all_genes if 'NDUFB7' in str(g).upper()]
    print(f"\nNDUFB7匹配 ({len(ndufb7_matches)}个):")
    for m in ndufb7_matches[:5]:
        print(f"  {m}")
    
    if not ndufb7_matches:
        # 尝试Ensembl ID
        ensg_matches = [g for g in all_genes if '99795' in str(g)]
        print(f"\nENSG99795匹配 ({len(ensg_matches)}个):")
        for m in ensg_matches[:5]:
            print(f"  {m}")
        target_genes = ensg_matches
    else:
        target_genes = ndufb7_matches
    
    if target_genes:
        target = target_genes[0]
        print(f"\n使用目标基因: {target}")
        
        results = []
        for fp in h5ad_files[:6]:
            try:
                adata = sc.read_h5ad(fp)
                
                # 修复空间坐标
                if 'X_spatial' in adata.obsm and 'spatial' not in adata.obsm:
                    adata.obsm['spatial'] = adata.obsm['X_spatial']
                
                if target not in adata.var_names:
                    print(f"[SKIP] {os.path.basename(fp)}: {target}不在基因列表")
                    continue
                
                sc.pp.normalize_total(adata, target_sum=1e4)
                sc.pp.log1p(adata)
                
                sq.gr.spatial_neighbors(adata)
                sq.gr.spatial_autocorr(adata, genes=[target], mode='moran')
                
                moran_i = float(adata.uns['moranI'][target]['I'])
                p_val = float(adata.uns['moranI'][target]['pval_norm'])
                
                # 推断区域
                region = "Unknown"
                for col in ['zone', 'region', 'sample', 'patient_region_id']:
                    if col in adata.obs.columns:
                        region = str(adata.obs[col].iloc[0])
                        break
                
                results.append({
                    'File': os.path.basename(fp),
                    'Region': region,
                    'N_Spots': adata.n_obs,
                    'NDUFB7_Mean': float(adata[:, target].X.mean()),
                    'Moran_I': moran_i,
                    'P_Value': p_val
                })
                print(f"[DONE] {os.path.basename(fp)}: MoranI={moran_i:.4f}, p={p_val:.2e}")
                
            except Exception as e:
                print(f"[ERROR] {os.path.basename(fp)}: {str(e)[:60]}")
        
        df = pd.DataFrame(results)
        df.to_csv(os.path.join(OUTDIR, "V107_Moran_I.csv"), index=False)
        print(f"\n[PASS] 保存: {OUTDIR}")
        print("\n结果汇总:")
        print(df.to_string(index=False))
    else:
        print("\n[FAIL] 未找到NDUFB7或ENSG00000099795")
PYEOF

echo ""
echo "========================================"
echo "V107 完成"
echo "========================================"
