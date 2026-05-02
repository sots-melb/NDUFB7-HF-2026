import os
import pandas as pd
import scipy.io
import numpy as np

project_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
results = []

print("▶ 开始全局巡检 GSE214611 (10x Matrix) 目录...")

# 1. 动态寻找所有包含 matrix.mtx 文件的目录
matrix_dirs = set()
for dp, dn, fns in os.walk(project_dir):
    if 'GSE214611' in dp and any('matrix.mtx' in f for f in fns):
        matrix_dirs.add(dp)

if not matrix_dirs:
    print("❌ 未找到 GSE214611 的 matrix.mtx 文件。请确保原始数据已解压。")
    exit(1)

for mdir in list(matrix_dirs):
    sample_name = os.path.basename(mdir)
    print(f"\n-> 正在解析样本: {sample_name}")

    # 2. 寻找特征文件 (处理可能的各种命名: features.tsv, genes.tsv, 带不带.gz)
    features_file = None
    for f in os.listdir(mdir):
        if 'features.tsv' in f or 'genes.tsv' in f:
            features_file = os.path.join(mdir, f)
            break

    if not features_file:
        print("   ❌ 找不到 features.tsv / genes.tsv, 跳过.")
        continue

    # 3. 无视格式，暴力检索 NDUFB7 所在行索引
    try:
        # header=None, on_bad_lines='skip' 防止格式破碎
        df_feat = pd.read_csv(features_file, sep='\t', header=None, on_bad_lines='skip')
        
        gene_idx = -1
        gene_name_found = ""
        # 遍历所有列，找 NDUFB7 (不区分大小写，兼容人类与小鼠)
        for col in df_feat.columns:
            if df_feat[col].dtype == object:
                matches = df_feat[df_feat[col].str.upper() == 'NDUFB7']
                if not matches.empty:
                    gene_idx = matches.index[0]
                    gene_name_found = matches.iloc[0][col]
                    break

        if gene_idx == -1:
            print("   ❌ 未在基因列表中找到 NDUFB7.")
            continue
            
        print(f"   ✅ 命中目标基因 [{gene_name_found}] (索引: {gene_idx}). 正在像手术刀一样切取矩阵...")

        # 4. 加载稀疏矩阵并直接切片
        mat_file = [f for f in os.listdir(mdir) if 'matrix.mtx' in f][0]
        mat_path = os.path.join(mdir, mat_file)
        
        # 读入为 CSR 格式加速行切片
        mat = scipy.io.mmread(mat_path).tocsr()
        
        # 10x 矩阵通常是 (genes x cells)，但也可能是反的
        if mat.shape[0] == len(df_feat):
            expr_array = mat[gene_idx, :].toarray()[0]
        elif mat.shape[1] == len(df_feat):
            expr_array = mat[:, gene_idx].toarray()[:, 0]
        else:
            print(f"   ❌ 矩阵维度 {mat.shape} 与基因数 {len(df_feat)} 不匹配!")
            continue

        # 5. 统计核心数据
        n_cells = len(expr_array)
        mean_expr = np.mean(expr_array)
        pct_zero = np.sum(expr_array == 0) / n_cells * 100

        print(f"   📊 细胞数: {n_cells} | NDUFB7均值: {mean_expr:.4f} | 零值率: {pct_zero:.2f}%")
        
        results.append({
            'Sample': sample_name,
            'Cells': n_cells,
            'Mean_NDUFB7': mean_expr,
            'Pct_Zero': pct_zero
        })

    except Exception as e:
        print(f"   ❌ 底层解析失败: {e}")

# 6. 保存最终结果
if results:
    df_res = pd.DataFrame(results)
    # 按样本名进行一个简单的时序排序尝试
    df_res = df_res.sort_values(by='Sample')
    
    out_csv = os.path.join(project_dir, "03_results/02_tables/GSE214611_Acute_Timeline_Rescued.csv")
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df_res.to_csv(out_csv, index=False)
    print(f"\n🎉 抢救大获成功！急/慢性时间轴已打通并保存至: {out_csv}")
    print("您可以直接查看结果：cat " + out_csv)
else:
    print("\n⚠️ 抢救未获取到任何数据。")
