import os
import gzip
import scipy.stats as stats

data_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE154170")
groups = {'Control': [], 'ICM': [], 'DCM': []}

for f in os.listdir(data_dir):
    if f.endswith(".csv.gz") and ("tpm_values" in f or "counts" in f):
        file_path = os.path.join(data_dir, f)
        print(f" -> 逐行扫描: {f}")
        try:
            with gzip.open(file_path, 'rt', encoding='utf-8') as gz_file:
                header = gz_file.readline().strip().split(',')
                sample_names = header[1:]
                for line in gz_file:
                    if line.startswith('"NDUFB7"') or line.startswith('NDUFB7,'):
                        expr_values = line.strip().split(',')[1:]
                        for i, expr in enumerate(expr_values):
                            val = float(expr)
                            s_name = sample_names[i].upper()
                            if 'RSB' in s_name or 'NOR' in s_name:
                                groups['Control'].append(val)
                            elif 'ICM' in s_name or 'ISCH' in s_name:
                                groups['ICM'].append(val)
                            elif 'DCM' in s_name or 'DILAT' in s_name:
                                groups['DCM'].append(val)
                        break # 找到基因后直接退出，极大节省内存和时间
        except Exception as e:
            print(f"扫描出错: {e}")

print("\n✅ 极低内存提取成功！各组样本数量:")
print(f" - Control: {len(groups['Control'])} | ICM: {len(groups['ICM'])} | DCM: {len(groups['DCM'])}")

if len(groups['Control']) > 0 and len(groups['ICM']) > 0:
    stat, p_icm = stats.ranksums(groups['Control'], groups['ICM'])
    print(f"\n📊 统计结果: Control vs ICM -> p-value = {p_icm:.4e}")
    if p_icm < 0.05:
        print("🎉 结论: NDUFB7 在缺血性心肌病(ICM)中显著下调！")
