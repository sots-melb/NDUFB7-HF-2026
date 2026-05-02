#!/bin/bash
# V104: F3表达提取 + A3 Moran's I重建
# 基于V103诊断确认的修复路径

PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"
DOWNLOADS="$HOME/Downloads"

echo "========================================"
echo "V104: F3 + A3 最终修复"
echo "========================================"

# ========================================
# PART 1: F3 GSE243655 featureCounts提取
# ========================================
echo ""
echo ">>> [F3] 从featureCounts提取NDUFB7"

EXTRACT_DIR="$PROJECT/01_data/01_raw_geo/GSE243655"
mkdir -p "$EXTRACT_DIR"

# 检查是否已解压
FC_FILES=$(ls "$EXTRACT_DIR"/*featureCounts.txt.gz 2>/dev/null | wc -l)

if [ "$FC_FILES" -lt 8 ]; then
    echo "  [EXTRACT] 解压RAW.tar到 $EXTRACT_DIR"
    tar -xf "$DOWNLOADS/GSE243655_RAW(1).tar" -C "$EXTRACT_DIR"
fi

echo "  featureCounts文件:"
ls "$EXTRACT_DIR"/*featureCounts.txt.gz | while read f; do
    echo "    $(basename "$f") ($(stat -c%s "$f" | numfmt --to=iec))"
done

# 提取NDUFB7表达值
echo ""
echo "  [EXTRACT] 逐样本提取NDUFB7..."

python3 << 'PYEOF'
import gzip
import pandas as pd
import os
from scipy import stats

EXTRACT_DIR = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE243655")

# 分组（已知）
groups = {
    'GSM7792365': 'DMSO', 'GSM7792366': 'DMSO', 
    'GSM7792367': 'DMSO', 'GSM7792368': 'DMSO',
    'GSM7792369': 'Fer-1', 'GSM7792370': 'Fer-1',
    'GSM7792371': 'Fer-1', 'GSM7792372': 'Fer-1'
}

results = []
for gsm, grp in groups.items():
    # 找文件（可能有不同命名前缀）
    import glob
    pattern = os.path.join(EXTRACT_DIR, f"{gsm}*featureCounts.txt.gz")
    files = glob.glob(pattern)
    
    if not files:
        print(f"[SKIP] {gsm} 文件未找到")
        continue
    
    fp = files[0]
    # featureCounts格式: Geneid Chr Start End Strand Length sample_counts
    df = pd.read_csv(fp, sep='\t', skiprows=1)  # skip header comment
    # 找NDUFB7行
    ndufb7_row = df[df.iloc[:,0].str.contains('NDUFB7', case=False, na=False)]
    
    if ndufb7_row.empty:
        # 尝试部分匹配
        ndufb7_row = df[df.iloc[:,0].str.contains('NDUF', case=False, na=False)]
        if not ndufb7_row.empty:
            print(f"  {gsm}: 找到NDUF家族基因: {ndufb7_row.iloc[0,0]}")
    
    if not ndufb7_row.empty:
        # 最后一列是counts
        counts = ndufb7_row.iloc[0, -1]
        gene_name = ndufb7_row.iloc[0, 0]
        results.append({
            'GSM': gsm, 'Group': grp, 'Gene': gene_name, 
            'Counts': int(counts), 'LogCounts': float(counts) + 1
        })
        print(f"  [OK] {gsm} ({grp}): {gene_name} = {counts}")
    else:
        print(f"  [FAIL] {gsm}: NDUFB7未找到")

# 统计
if len(results) == 8:
    df_res = pd.DataFrame(results)
    dmso = df_res[df_res['Group']=='DMSO']['Counts'].values
    fer1 = df_res[df_res['Group']=='Fer-1']['Counts'].values
    
    # 配对t-test（同患者配对）
    t_stat, p_val = stats.ttest_rel(dmso, fer1)
    delta = fer1.mean() - dmso.mean()
    pct = delta / dmso.mean() * 100
    
    print(f"\n=== F3 配对t-test结果 ===")
    print(f"DMSO mean counts: {dmso.mean():.1f}")
    print(f"Fer-1 mean counts: {fer1.mean():.1f}")
    print(f"Δ = {delta:.1f} ({pct:+.1f}%)")
    print(f"t = {t_stat:.3f}")
    print(f"p = {p_val:.4f}")
    
    if p_val < 0.05:
        if fer1.mean() > dmso.mean():
            print(f"\n[PASS] Fer-1显著上调NDUFB7！GRADE A达成")
        else:
            print(f"\n[PASS] 显著但方向需确认")
    else:
        print(f"\n[WARN] 不显著 (p={p_val:.3f})")
    
    # 保存
    outdir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/03_results/V104_F3_Fer1")
    os.makedirs(outdir, exist_ok=True)
    df_res.to_csv(os.path.join(outdir, "V104_Fer1_NDUFB7_featureCounts.csv"), index=False)
    
    # 配对结果
    paired = pd.DataFrame({
        'Patient': [1,2,3,4],
        'DMSO': dmso,
        'Fer1': fer1,
        'Delta': fer1 - dmso
    })
    paired.to_csv(os.path.join(outdir, "V104_Fer1_paired.csv"), index=False)
    print(f"\n[DONE] 保存: {outdir}")
else:
    print(f"\n[FAIL] 仅提取到 {len(results)}/8 样本")
PYEOF

# ========================================
# PART 2: A3 Moran's I重建（X_spatial修复）
# ========================================
echo ""
echo ">>> [A3] Moran's I重建（X_spatial键修复）"

python3 << 'PYEOF'
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import os
from datetime import datetime

PROJECT_DIR = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "03_results/V104_A3_Moran_I")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 扫描所有Visium h5ad
h5ad_files = []
for root, dirs, files in os.walk(os.path.join(PROJECT_DIR, "01_data")):
    for f in files:
        if f.endswith('.h5ad') and 'Visium' in f:
            h5ad_files.append(os.path.join(root, f))

print(f"[{datetime.now().strftime('%H:%M:%S')}] === V104 A3: Moran's I重建 ===")
print(f"找到 {len(h5ad_files)} 个Visium文件")

results = []

for fp in h5ad_files[:6]:  # 最多6个
    try:
        adata = sc.read_h5ad(fp)
        print(f"\n处理: {os.path.basename(fp)}")
        print(f"  维度: {adata.n_obs} spots × {adata.n_vars} genes")
        
        # 修复：将X_spatial复制为spatial（Squidpy默认找spatial）
        if 'X_spatial' in adata.obsm and 'spatial' not in adata.obsm:
            adata.obsm['spatial'] = adata.obsm['X_spatial']
            print(f"  [FIX] 复制 X_spatial -> spatial")
        
        if 'spatial' not in adata.obsm:
            print(f"  [SKIP] 无空间坐标")
            continue
        
        # 质控+标准化
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # 检查NDUFB7
        if 'NDUFB7' not in adata.var_names:
            print(f"  [WARN] NDUFB7不在基因列表")
            continue
        
        # 计算空间邻域图
        sq.gr.spatial_neighbors(adata)
        
        # 计算Moran's I
        sq.gr.spatial_autocorr(adata, genes=['NDUFB7'], mode='moran')
        
        moran_i = float(adata.uns['moranI']['NDUFB7']['I'])
        p_val = float(adata.uns['moranI']['NDUFB7']['pval_norm'])
        
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
            'NDUFB7_Mean': float(adata[:, 'NDUFB7'].X.mean()),
            'NDUFB7_Pct': float((adata[:, 'NDUFB7'].X > 0).mean() * 100),
            'Moran_I': moran_i,
            'P_Value': p_val,
            'Status': 'Success'
        })
        
        print(f"  [DONE] Moran's I = {moran_i:.4f}, p = {p_val:.2e}")
        
    except Exception as e:
        print(f"  [ERROR] {str(e)[:80]}")
        results.append({
            'File': os.path.basename(fp),
            'Region': 'Error',
            'N_Spots': 0,
            'NDUFB7_Mean': np.nan,
            'NDUFB7_Pct': np.nan,
            'Moran_I': np.nan,
            'P_Value': np.nan,
            'Status': f'Error: {str(e)[:40]}'
        })

# 保存
df = pd.DataFrame(results)
csv_path = os.path.join(OUTPUT_DIR, f"V104_Moran_I_{datetime.now().strftime('%Y%m%d')}.csv")
df.to_csv(csv_path, index=False)

print(f"\n[{datetime.now().strftime('%H:%M:%S')}] === 结果汇总 ===")
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
echo "V104 完成"
echo "========================================"
echo "F3结果: 03_results/V104_F3_Fer1/"
echo "A3结果: 03_results/V104_A3_Moran_I/"
