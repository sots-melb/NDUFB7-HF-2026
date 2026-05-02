import gzip
import pandas as pd
import os

in_file = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz")
out_file = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE59867_Clinical_Raw.csv")

print(f"正在读取 {in_file}...")
try:
    lines = []
    with gzip.open(in_file, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_'):
                lines.append(line.strip().split('\t'))
                
    # 转置并清洗数据
    df = pd.DataFrame(lines).set_index(0).T
    df.columns = [col.replace('!Sample_', '') for col in df.columns]
    df.to_csv(out_file, index=False)
    print(f"✅ 临床数据已暴力提取并转置成功！共 {len(df)} 个样本。保存至 {out_file}")
except Exception as e:
    print(f"❌ 提取失败: {e}")
