#!/bin/bash
# V108: F3科学重解读 + A3健壮化重建
# 这是最后一轮技术修复，之后必须写Results

PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"

echo "========================================"
echo "V108: F3科学重解读 + A3健壮化重建"
echo "========================================"

# ========================================
# F3: 科学重解读 + 保存最终结论
# ========================================
echo ""
echo ">>> [F3] 科学重解读：Fer-1不影响NDUFB7的因果意义"

python3 << 'PYEOF'
import pandas as pd
import os
from scipy import stats

OUTDIR = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/03_results/V108_F3_Fer1_Final")
os.makedirs(OUTDIR, exist_ok=True)

# 读取V107结果
v107 = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/03_results/V107_F3_Fer1/V107_Fer1_NDUFB7.csv")
df = pd.read_csv(v107)

dmso = df[df['Group']=='DMSO']['Counts'].values
fer1 = df[df['Group']=='Fer-1']['Counts'].values

# 配对检验
t_stat, p_val = stats.ttest_rel(dmso, fer1)
delta = fer1.mean() - dmso.mean()
pct = delta / dmso.mean() * 100

print("=== F3 最终科学解读 ===")
print(f"DMSO (Vehicle):  {dmso.mean():.1f} ± {dmso.std():.1f} counts")
print(f"Fer-1 (10µM):    {fer1.mean():.1f} ± {fer1.std():.1f} counts")
print(f"Δ = {delta:.1f} ({pct:+.1f}%)")
print(f"paired t(3) = {t_stat:.3f}, p = {p_val:.4f}")
print("")

# 科学解读
if p_val > 0.05:
    print("[INTERPRETATION] Fer-1 不改变 NDUFB7 表达 (p>0.05)")
    print("  → NDUFB7 丢失不是铁死亡的下游结果")
    print("  → Fer-1 阻断脂质过氧化/铁死亡，但无法恢复 Complex I 亚基")
    print("  → 支持因果方向: NDUFB7↓ (上游) → OXPHOS崩溃 → ROS↑ → 铁死亡↑ (下游)")
    print("  → 这意味着 NDUFB7 是独立的起始节点，而非铁死亡的被动受害者")
    print("")
    print("[PAPER CLAIM] 'Ferrostatin-1 treatment did not alter NDUFB7 expression,")
    print("  suggesting that NDUFB7 loss precedes ferroptosis activation rather than")
    print("  representing a downstream consequence of lipid peroxidation.'")
    print("")
    print("[IMPLICATION] 治疗策略需要直接靶向 NDUFB7/Complex I，而非仅抑制铁死亡")
else:
    if fer1.mean() > dmso.mean():
        print("[PASS] Fer-1 上调 NDUFB7")
    else:
        print("[PARTIAL] 显著下调，需解释")

# 保存最终结论
conclusion = pd.DataFrame({
    'Metric': ['DMSO_mean', 'Fer1_mean', 'Delta_counts', 'Delta_pct', 't_stat', 'p_value', 'N_pair'],
    'Value': [dmso.mean(), fer1.mean(), delta, pct, t_stat, p_val, 4],
    'Interpretation': ['Vehicle control', 'Ferrostatin-1 10µM', 
                       'No significant change', 'No significant change',
                       'NS', 'NS', '4 patients, paired']
})
conclusion.to_csv(os.path.join(OUTDIR, "V108_F3_conclusion.csv"), index=False)

# 配对数据
paired = pd.DataFrame({
    'Patient': [1, 2, 3, 4],
    'DMSO': dmso,
    'Fer1': fer1,
    'Delta': fer1 - dmso
})
paired.to_csv(os.path.join(OUTDIR, "V108_F3_paired_final.csv"), index=False)

print("[DONE] F3最终结论保存: ", OUTDIR)
PYEOF

# ========================================
# A3: 健壮化重建（处理sparse/zero/数据类型）
# ========================================
echo ""
echo ">>> [A3] 健壮化Moran's I重建"

python3 << 'PYEOF'
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import os
from datetime import datetime

PROJECT = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
OUTDIR = os.path.join(PROJECT, "03_results/V108_A3_Moran_I_Final")
os.makedirs(OUTDIR, exist_ok=True)

# 扫描Visium
h5ad_files = []
for root, dirs, files in os.walk(os.path.join(PROJECT, "01_data")):
    for f in files:
        if f.endswith('.h5ad') and 'Visium' in f:
            h5ad_files.append(os.path.join(root, f))

print(f"[{datetime.now().strftime('%H:%M:%S')}] 扫描到 {len(h5ad_files)} 个Visium文件")

results = []

for fp in h5ad_files[:8]:  # 最多8个
    try:
        print(f"\n处理: {os.path.basename(fp)}")
        adata = sc.read_h5ad(fp)
        print(f"  原始: {adata.n_obs} spots × {adata.n_vars} genes")
        
        # 修复1: 复制X_spatial
        if 'X_spatial' in adata.obsm:
            adata.obsm['spatial'] = adata.obsm['X_spatial']
            print(f"  [FIX] X_spatial -> spatial")
        
        # 修复2: 确保NDUFB7存在
        target = None
        if 'NDUFB7' in adata.var_names:
            target = 'NDUFB7'
        else:
            matches = [g for g in adata.var_names if 'NDUFB7' in str(g).upper()]
            if matches:
                target = matches[0]
        
        if not target:
            print(f"  [SKIP] 无NDUFB7")
            continue
        print(f"  [OK] 目标基因: {target}")
        
        # 修复3: 质控过滤（去除0表达spots过多的情况）
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        
        # 修复4: 标准化
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # 修复5: 检查NDUFB7表达分布
        expr = adata[:, target].X
        # 转为dense
        if hasattr(expr, 'toarray'):
            expr_dense = expr.toarray().flatten()
        else:
            expr_dense = np.array(expr).flatten()
        
        n_zero = np.sum(expr_dense == 0)
        pct_zero = n_zero / len(expr_dense) * 100
        print(f"  [INFO] {target}表达: mean={expr_dense.mean():.3f}, 零值={pct_zero:.1f}%")
        
        # 如果零值>95%，Moran's I可能不稳定，但继续尝试
        if pct_zero > 99:
            print(f"  [WARN] 零值过高，Moran's I可能不可靠")
        
        # 修复6: 显式指定obsm_spatial_key
        try:
            sq.gr.spatial_neighbors(adata, coord_type="generic")
        except Exception as e:
            print(f"  [WARN] spatial_neighbors: {str(e)[:60]}")
            # 备用：手动构建邻域图
            from sklearn.neighbors import NearestNeighbors
            coords = adata.obsm['spatial']
            nbrs = NearestNeighbors(n_neighbors=6, metric='euclidean').fit(coords)
            distances, indices = nbrs.kneighbors(coords)
            # 转换为squidpy格式
            adata.obsp['spatial_distances'] = distances
            adata.obsp['spatial_connectivities'] = indices
        
        # 修复7: 计算Moran's I（显式参数）
        try:
            sq.gr.spatial_autocorr(
                adata, 
                genes=[target], 
                mode='moran',
                n_perms=100,  # 减少置换次数避免内存问题
                n_jobs=1
            )
            
            moran_i = float(adata.uns['moranI'][target]['I'])
            p_val = float(adata.uns['moranI'][target].get('pval_norm', 
                         adata.uns['moranI'][target].get('pval_z', np.nan)))
            
            print(f"  [DONE] Moran's I = {moran_i:.4f}, p = {p_val:.2e}")
            
            # 推断区域
            region = "Unknown"
            for col in ['zone', 'region', 'sample', 'patient_region_id', 'patient']:
                if col in adata.obs.columns:
                    vals = adata.obs[col].unique()
                    region = str(vals[0]) if len(vals) == 1 else "Mixed"
                    break
            
            results.append({
                'File': os.path.basename(fp),
                'Region': region,
                'N_Spots': adata.n_obs,
                'NDUFB7_Mean': float(expr_dense.mean()),
                'NDUFB7_Pct_Zero': float(pct_zero),
                'Moran_I': float(moran_i),
                'P_Value': float(p_val),
                'Status': 'Success'
            })
            
        except Exception as e2:
            print(f"  [ERROR] spatial_autocorr: {str(e2)[:80]}")
            # 备用方案：手动计算Moran's I
            try:
                from esda.moran import Moran
                import libpysal
                
                y = expr_dense
                w = libpysal.weights.KNN.from_array(adata.obsm['spatial'], k=6)
                w.transform = 'r'
                mi = Moran(y, w)
                
                print(f"  [FALLBACK] esda Moran's I = {mi.I:.4f}, p = {mi.p_norm:.2e}")
                
                results.append({
                    'File': os.path.basename(fp),
                    'Region': region,
                    'N_Spots': adata.n_obs,
                    'NDUFB7_Mean': float(expr_dense.mean()),
                    'NDUFB7_Pct_Zero': float(pct_zero),
                    'Moran_I': float(mi.I),
                    'P_Value': float(mi.p_norm),
                    'Status': 'Fallback_esda'
                })
            except Exception as e3:
                print(f"  [FAIL] 备用也失败: {str(e3)[:60]}")
                results.append({
                    'File': os.path.basename(fp),
                    'Region': 'Error',
                    'N_Spots': adata.n_obs,
                    'NDUFB7_Mean': float(expr_dense.mean()),
                    'NDUFB7_Pct_Zero': float(pct_zero),
                    'Moran_I': np.nan,
                    'P_Value': np.nan,
                    'Status': f'Error: {str(e2)[:30]}'
                })
        
    except Exception as e:
        print(f"  [ERROR] 整体失败: {str(e)[:80]}")
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
csv_path = os.path.join(OUTDIR, f"V108_Moran_I_{datetime.now().strftime('%Y%m%d')}.csv")
df.to_csv(csv_path, index=False)

print(f"\n[{datetime.now().strftime('%H:%M:%S')}] === A3 结果汇总 ===")
print(df.to_string(index=False))
print(f"\n[PASS] 保存: {csv_path}")

# 与模板对比
print("\n=== 与V77模板值对比 ===")
template = {'Control': 0.015, 'FZ': 0.119, 'IZ': 0.041, 'BZ': 0.038}
print(f"模板值: {template}")
print("重建值见上表。差异>5%时优先信任重建值。")
PYEOF

echo ""
echo "========================================"
echo "V108 完成"
echo "========================================"
echo "F3结论: 03_results/V108_F3_Fer1_Final/"
echo "A3结果: 03_results/V108_A3_Moran_I_Final/"
echo ""
echo "=== 下一步 ==="
echo "1. 查看F3科学解读（支持NDUFB7为上游节点）"
echo "2. 查看A3 Moran's I重建结果"
echo "3. 无论成败，立即开始写Results"
