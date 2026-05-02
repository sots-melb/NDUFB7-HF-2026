import os
import pandas as pd
import scipy.io
import numpy as np

project_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE214611_FULL_EXTRACT")
results = []

print("\n▶ 开始在完整目录巡检 10x Matrix...")

matrix_dirs = set()
for dp, dn, fns in os.walk(project_dir):
    if any('matrix.mtx' in f for f in fns):
        matrix_dirs.add(dp)

for mdir in list(matrix_dirs):
    sample_name = os.path.basename(mdir)
    
    features_file = None
    for f in os.listdir(mdir):
        if 'features.tsv' in f or 'genes.tsv' in f:
            features_file = os.path.join(mdir, f)
            break

    if not features_file: continue

    try:
        df_feat = pd.read_csv(features_file, sep='\t', header=None, on_bad_lines='skip')
        gene_idx = -1
        for col in df_feat.columns:
            if df_feat[col].dtype == object:
                matches = df_feat[df_feat[col].str.upper() == 'NDUFB7']
                if not matches.empty:
                    gene_idx = matches.index[0]
                    break
        if gene_idx == -1: continue

        mat_file = [f for f in os.listdir(mdir) if 'matrix.mtx' in f][0]
        mat_path = os.path.join(mdir, mat_file)
        mat = scipy.io.mmread(mat_path).tocsr()
        
        if mat.shape[0] == len(df_feat): expr_array = mat[gene_idx, :].toarray()[0]
        elif mat.shape[1] == len(df_feat): expr_array = mat[:, gene_idx].toarray()[:, 0]
        else: continue

        n_cells = len(expr_array)
        mean_expr = np.mean(expr_array)
        
        # 仅打印找到了有效表达的样本
        if mean_expr > 0:
            print(f"   ✅ [样本: {sample_name:25s}] 细胞数: {n_cells:<6d} | 均值: {mean_expr:.4f}")
            results.append({'Sample': sample_name, 'Cells': n_cells, 'Mean_NDUFB7': mean_expr})

    except Exception as e:
        pass

if results:
    df_res = pd.DataFrame(results).sort_values(by='Sample')
    out_csv = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE214611_Full_Acute_Timeline.csv")
    df_res.to_csv(out_csv, index=False)
    print(f"\n🎉 终极完整时间轴提取完毕！已保存至: {out_csv}")
