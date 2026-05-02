import scanpy as sc
import squidpy as sq
import os

project_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
h5ad_paths = [os.path.join(dp, f) for dp, dn, filenames in os.walk(project_dir) for f in filenames if f.endswith('.h5ad') and 'Kuppe' in dp]

if h5ad_paths:
    # backed='r' 模式：只读磁盘不吃内存
    try:
        adata = sc.read_h5ad(h5ad_paths[0], backed='r')
    except:
        adata = sc.read_h5ad(h5ad_paths[0])
        
    spatial_key = 'spatial'
    if 'spatial' not in adata.obsm.keys():
        possible_keys = [k for k in adata.obsm.keys() if 'spatial' in k.lower()]
        if possible_keys: spatial_key = possible_keys[0]
            
    sq.gr.spatial_neighbors(adata, spatial_key=spatial_key)
    target_gene = next((g for g in adata.var_names if g.upper() == 'NDUFB7'), None)
    
    if target_gene:
        print(f"✅ 找到 {target_gene}，开始计算 Moran's I...")
        if adata.isbacked:
            adata = adata[:, [target_gene]].to_memory()
            sq.gr.spatial_neighbors(adata, spatial_key=spatial_key)
            
        sq.gr.spatial_autocorr(adata, mode="moran", genes=[target_gene])
        
        df_res = adata.uns['moranI']
        moran_i_val = df_res.loc[target_gene, 'I']
        p_val = df_res.loc[target_gene, 'pval_sim']
        
        print("\n==========================================")
        print(f"🌟 NDUFB7 真实 Moran's I: {moran_i_val:.4f}")
        print(f"🌟 显著性 p-value: {p_val:.4e}")
        print("==========================================")
