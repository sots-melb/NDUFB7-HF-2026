#!/usr/bin/env python3
import pandas as pd
import sys

gtex_dir = "/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/gtex/GTEx_Analysis_v11_eQTL"
parquet_file = f"{gtex_dir}/Heart_Left_Ventricle.v11.eQTLs.signif_pairs.parquet"

print("="*70)
print("GTEx v11 Heart LV: NDUFB7 parquet提取")
print("="*70)

try:
    # 读取全部数据（内存允许的话）
    df = pd.read_parquet(parquet_file)
    print(f"总signif_pairs: {len(df)}")
    print(f"列名: {list(df.columns)}")
    
    # 提取NDUFB7
    ndufb7 = df[df['gene_id'].str.contains('ENSG00000099795', na=False)]
    print(f"\nNDUFB7 signif_pairs: {len(ndufb7)}")
    
    if len(ndufb7) > 0:
        print("\n前5行:")
        print(ndufb7.head())
        
        # 保存
        out_file = "/home/y411869/Projects/NDUFB7_HF_2026_04_20/03_results/gtex_ndufb7_all_pairs.txt"
        ndufb7.to_csv(out_file, sep='\t', index=False)
        print(f"\n✅ 已保存: {out_file}")
        
        # 统计
        print(f"\n【统计】")
        print(f"  唯一SNP数: {ndufb7['variant_id'].nunique()}")
        print(f"  slope范围: {ndufb7['slope'].min():.4f} ~ {ndufb7['slope'].max():.4f}")
        print(f"  pval_nominal范围: {ndufb7['pval_nominal'].min():.2e} ~ {ndufb7['pval_nominal'].max():.2e}")
        
except Exception as e:
    print(f"❌ 错误: {e}")
    print("可能需要安装: pip install pyarrow fastparquet")
    sys.exit(1)
