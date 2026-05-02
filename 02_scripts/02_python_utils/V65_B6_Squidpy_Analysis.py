import scanpy as sc
import squidpy as sq
import os

print("启动 Squidpy 空间邻域分析...")
project_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
# 动态寻找 Kuppe Visium .h5ad 文件 (S01 资产)
h5ad_paths = [os.path.join(dp, f) for dp, dn, filenames in os.walk(project_dir) for f in filenames if f.endswith('.h5ad') and 'Kuppe' in dp]

if h5ad_paths:
    adata = sc.read_h5ad(h5ad_paths[0])
    print(f"✅ 成功加载 Visium 数据: {h5ad_paths[0]}")
    
    # 构建空间图并计算邻域富集
    sq.gr.spatial_neighbors(adata)
    # 假设细胞类型注释在 'cell_type' 列，如果不存在则跳过实际计算以防报错
    if 'cell_type' in adata.obs.columns:
        sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
        sq.gr.co_occurrence(adata, cluster_key="cell_type")
        print("✅ 空间邻域 (nhood_enrichment) 与共现 (co_occurrence) 计算完成。")
    else:
        print("⚠️ 未找到 cell_type 列，等待 AddModuleScore / 去卷积完成后再运行后续绘图。")
else:
    print("⚠️ 未找到 Kuppe Visium 的 .h5ad 文件，需确认数据位置。")
