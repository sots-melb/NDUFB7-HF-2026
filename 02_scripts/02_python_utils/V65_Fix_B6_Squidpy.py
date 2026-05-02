import scanpy as sc
import squidpy as sq
import os

print("启动 Squidpy 空间邻域分析 (自适应坐标键名)...")
project_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
h5ad_paths = [os.path.join(dp, f) for dp, dn, filenames in os.walk(project_dir) for f in filenames if f.endswith('.h5ad') and 'Kuppe' in dp]

if h5ad_paths:
    adata = sc.read_h5ad(h5ad_paths[0])
    print(f"✅ 成功加载 Visium 数据: {h5ad_paths[0]}")
    
    # 动态探查空间坐标键名
    spatial_key = 'spatial'
    if 'spatial' not in adata.obsm.keys():
        possible_keys = [k for k in adata.obsm.keys() if 'spatial' in k.lower()]
        if possible_keys:
            spatial_key = possible_keys[0]
            print(f"🔍 动态识别到空间坐标 key 为: {spatial_key}")
        else:
            print(f"❌ 未在 adata.obsm 中找到包含 spatial 的键: {list(adata.obsm.keys())}")
            exit(1)
            
    # 构建空间图
    sq.gr.spatial_neighbors(adata, spatial_key=spatial_key)
    print("✅ 空间网络构建成功，Squidpy 环境及数据完全畅通！")
else:
    print("⚠️ 未找到 Kuppe Visium 的 .h5ad 文件。")
