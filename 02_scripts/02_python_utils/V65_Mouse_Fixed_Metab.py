import os
import pandas as pd
import scipy.io
import numpy as np
import scanpy as sc
import scipy.stats as stats
import warnings
warnings.filterwarnings('ignore')

project_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE214611_FULL_EXTRACT")

# 定位 snd3_2 样本 (心梗后 3 天，处于机制活跃期)
target_sample = next((dp for dp, dn, fns in os.walk(project_dir) if 'snd3_2' in dp and any('matrix.mtx' in f for f in fns)), None)

if target_sample:
    mat_file = [f for f in os.listdir(target_sample) if 'matrix.mtx' in f][0]
    feat_file = [f for f in os.listdir(target_sample) if 'features.tsv' in f or 'genes.tsv' in f][0]
    
    df_feat = pd.read_csv(os.path.join(target_sample, feat_file), sep='\t', header=None, on_bad_lines='skip')
    
    gene_col = next((col for col in df_feat.columns if df_feat[col].dtype == object and 'Ndufb7' in df_feat[col].values), 1)
    gene_names = df_feat[gene_col].values
    
    mat = scipy.io.mmread(os.path.join(target_sample, mat_file)).tocsr()
    if mat.shape[0] == len(gene_names): mat = mat.T
        
    adata = sc.AnnData(X=mat)
    adata.var_names = gene_names
    adata.var_names_make_unique()
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # --- 1. 绝对提纯心肌细胞 ---
    cm_genes = [g for g in ['Tnnt2', 'Myh6', 'Myh7'] if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=cm_genes, score_name='CM_Score')
    adata_cm = adata[adata.obs['CM_Score'] > 0.5].copy()
    
    print(f"✅ 成功提纯小鼠心肌细胞: {adata_cm.n_obs} 个")
    
    # --- 2. 机制打分 ---
    oxphos = [g for g in ['Atp5f1a', 'Cox4i1', 'Cyc1', 'Ndufa4', 'Sdha'] if g in adata_cm.var_names]
    ferro_def = [g for g in ['Slc7a11', 'Gpx4', 'Fth1', 'Ftl1', 'Nfe2l2', 'Gss'] if g in adata_cm.var_names]
    ferro_drv = [g for g in ['Acsl4', 'Ptgs2', 'Alox15', 'Tfrc'] if g in adata_cm.var_names]
    
    sc.tl.score_genes(adata_cm, gene_list=oxphos, score_name='OXPHOS')
    sc.tl.score_genes(adata_cm, gene_list=ferro_def, score_name='Ferro_Def')
    sc.tl.score_genes(adata_cm, gene_list=ferro_drv, score_name='Ferro_Drv')
    
    # --- 3. 统计相关性 (修复报错核心) ---
    if 'Ndufb7' in adata_cm.var_names:
        x_data = adata_cm[:, 'Ndufb7'].X
        # 修复点：强制解压稀疏矩阵
        expr = x_data.toarray().flatten() if hasattr(x_data, 'toarray') else np.array(x_data).flatten()
        
        r_ox, p_ox = stats.pearsonr(expr, adata_cm.obs['OXPHOS'])
        r_def, p_def = stats.pearsonr(expr, adata_cm.obs['Ferro_Def'])
        r_drv, p_drv = stats.pearsonr(expr, adata_cm.obs['Ferro_Drv'])
        
        print("\n=======================================================")
        print("🌟 机制闭环 3.0：小鼠单核纯 CM 真实相关性 🌟")
        print(f" 1. Ndufb7 vs OXPHOS 产能:       r = {r_ox:.4f} (p={p_ox:.2e})")
        print(f" 2. Ndufb7 vs 铁死亡防线(Gpx4):  r = {r_def:.4f} (p={p_def:.2e})")
        print(f" 3. Ndufb7 vs 铁死亡驱动(Acsl4): r = {r_drv:.4f} (p={p_drv:.2e})")
        print("=======================================================")
    else:
        print("❌ 未找到 Ndufb7")
else:
    print("❌ 未找到 snd3_2 样本")
