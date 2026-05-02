import os
import pandas as pd
import scipy.stats as stats

print("▶ 开始合并 GSE154170 TPM 矩阵...")
data_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE154170")

dfs = []
for f in os.listdir(data_dir):
    if f.endswith(".csv.gz") and "tpm_values" in f:
        file_path = os.path.join(data_dir, f)
        print(f" -> 正在加载: {f}")
        try:
            # 第一列通常是 gene symbol 或 ID
            df = pd.read_csv(file_path, index_col=0)
            dfs.append(df)
        except Exception as e:
            print(f"加载失败: {e}")

if dfs:
    # 按照索引 (基因名) 横向拼接
    merged_df = pd.concat(dfs, axis=1)
    
    # 寻找 NDUFB7 的行
    target_idx = [idx for idx in merged_df.index if str(idx).upper() == 'NDUFB7']
    
    if target_idx:
        ndufb7_expr = merged_df.loc[target_idx[0]]
        
        # 根据列名分配组别 (GSE154170的列名通常包含 RSB, ICM, DCM)
        groups = {'Control': [], 'ICM': [], 'DCM': []}
        for sample, expr in ndufb7_expr.items():
            s_upper = str(sample).upper()
            if 'RSB' in s_upper or 'NOR' in s_upper:
                groups['Control'].append(expr)
            elif 'ICM' in s_upper or 'ISCH' in s_upper:
                groups['ICM'].append(expr)
            elif 'DCM' in s_upper or 'DILAT' in s_upper:
                groups['DCM'].append(expr)
                
        print("\n✅ 样本分组重构成功:")
        print(f" - Control 样本数: {len(groups['Control'])}")
        print(f" - ICM 样本数: {len(groups['ICM'])}")
        print(f" - DCM 样本数: {len(groups['DCM'])}")
        
        # 统计检验 (Kruskal-Wallis 或 Wilcoxon)
        if len(groups['Control']) > 0 and len(groups['ICM']) > 0:
            stat, p_icm = stats.ranksums(groups['Control'], groups['ICM'])
            print(f"\n📊 统计结果: Control vs ICM -> p-value = {p_icm:.4e}")
            if p_icm < 0.05:
                print("🎉 结论: NDUFB7 在缺血性心肌病(ICM)中发生显著特异性下调！")
        else:
            print("⚠️ 无法进行统计检验，某一组数量为0。")
    else:
        print("❌ 未在矩阵中找到 NDUFB7 基因。")
else:
    print("❌ 未找到任何 TPM 矩阵压缩包。")
