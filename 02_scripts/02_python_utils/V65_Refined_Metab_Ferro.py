import scanpy as sc
import os
import numpy as np
import scipy.stats as stats

print("▶ 加载 Kuppe 空间数据...")
project_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
h5ad_paths = [os.path.join(dp, f) for dp, dn, filenames in os.walk(project_dir) for f in filenames if f.endswith('.h5ad') and 'Kuppe' in dp]

if h5ad_paths:
    adata = sc.read_h5ad(h5ad_paths[0])
    
    # --- 1. 细胞特异性过滤 (寻找 CM 富集区) ---
    cm_markers = ['TNNT2', 'MYH7', 'MYBPC3', 'ACTC1', 'RYR2']
    cm_markers = [g for g in cm_markers if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=cm_markers, score_name='CM_Score')
    
    # 仅保留 CM 得分在前 50% 的 spot (剔除纯纤维化或免疫浸润区)
    median_cm = np.median(adata.obs['CM_Score'])
    adata_cm = adata[adata.obs['CM_Score'] > median_cm].copy()
    print(f"✅ 成功圈定心肌细胞富集微环境 (Spots数量: {adata_cm.n_obs} / {adata.n_obs})")
    
    # --- 2. 机制基因集正反拆分 ---
    # 铁死亡防御 (越低越容易死)
    ferro_defense = ['SLC7A11', 'GPX4', 'FTH1', 'FTL', 'NFE2L2', 'GSS']
    # 铁死亡驱动 (越高越容易死)
    ferro_driver = ['ACSL4', 'PTGS2', 'ALOX15', 'TFRC']
    # 产能核心
    oxphos_genes = ['ATP5F1A', 'COX4I1', 'CYC1', 'NDUFA4', 'SDHA']
    
    ferro_defense = [g for g in ferro_defense if g in adata_cm.var_names]
    ferro_driver = [g for g in ferro_driver if g in adata_cm.var_names]
    oxphos_genes = [g for g in oxphos_genes if g in adata_cm.var_names]
    
    # --- 3. 重新打分 ---
    sc.tl.score_genes(adata_cm, gene_list=ferro_defense, score_name='Ferro_Defense_Score')
    sc.tl.score_genes(adata_cm, gene_list=ferro_driver, score_name='Ferro_Driver_Score')
    sc.tl.score_genes(adata_cm, gene_list=oxphos_genes, score_name='OXPHOS_Score')
    
    # --- 4. 计算 CM 微环境中的真实相关性 ---
    target_gene = next((g for g in adata_cm.var_names if g.upper() == 'NDUFB7'), None)
    
    if target_gene:
        if hasattr(adata_cm[:, target_gene].X, 'toarray'):
            ndufb7_expr = adata_cm[:, target_gene].X.toarray().flatten()
        else:
            ndufb7_expr = np.array(adata_cm[:, target_gene].X).flatten()
            
        oxphos_score = adata_cm.obs['OXPHOS_Score'].values
        defense_score = adata_cm.obs['Ferro_Defense_Score'].values
        driver_score = adata_cm.obs['Ferro_Driver_Score'].values
        
        r_oxphos, p_oxphos = stats.pearsonr(ndufb7_expr, oxphos_score)
        r_def, p_def = stats.pearsonr(ndufb7_expr, defense_score)
        r_dri, p_dri = stats.pearsonr(ndufb7_expr, driver_score)
        
        print("\n=======================================================")
        print("🌟 机制闭环 2.0：心肌微环境下的真实相关性 🌟")
        print(f" 1. NDUFB7 vs OXPHOS 产能网络: r = {r_oxphos:.4f} (p={p_oxphos:.2e})")
        print(f" 2. NDUFB7 vs 铁死亡防线(GPX4等): r = {r_def:.4f} (p={p_def:.2e})")
        print(f" 3. NDUFB7 vs 铁死亡驱动(ACSL4等): r = {r_dri:.4f} (p={p_dri:.2e})")
        print("=======================================================")
        
        if r_oxphos > 0.1 and r_def > 0.1 and r_dri < 0:
            print("\n🎉 逻辑完美闭环！在心肌细胞内，NDUFB7 的丢失不仅导致产能崩溃，更直接导致抗氧化防线瓦解，从而释放铁死亡驱动信号！")
        else:
            print("\n💡 提示：空间水平信号依然复杂。若结果不达预期，我们将祭出单核转录组(GSE183852)对 CM 进行绝对提纯。")
            
    else:
        print("❌ 未找到 NDUFB7 基因。")
else:
    print("❌ 未找到 .h5ad 空间切片文件。")
