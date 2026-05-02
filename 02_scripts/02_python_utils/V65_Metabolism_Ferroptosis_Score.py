import scanpy as sc
import os
import numpy as np
import scipy.stats as stats

print("▶ 加载 Kuppe 空间数据进行代谢评分...")
project_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
h5ad_paths = [os.path.join(dp, f) for dp, dn, filenames in os.walk(project_dir) for f in filenames if f.endswith('.h5ad') and 'Kuppe' in dp]

if h5ad_paths:
    # 完整加载入内存进行计算 (空间数据不大)
    adata = sc.read_h5ad(h5ad_paths[0])
    
    # --- 构建机制基因集 ---
    # 1. 核心铁死亡防线 (V64 确立)
    ferro_defense = ['FTL', 'SAT1', 'NFE2L2', 'SLC7A11', 'ACSL4']
    # 2. 氧化磷酸化 (OXPHOS) 核心产能亚基
    oxphos_genes = ['ATP5F1A', 'COX4I1', 'CYC1', 'NDUFA4', 'SDHA']
    # 3. 糖酵解 (Glycolysis) 核心酶 (测试是否有代偿)
    glyco_genes = ['HK1', 'HK2', 'PFKM', 'ALDOA', 'ENO1', 'PKM', 'LDHA']
    
    # 确保基因在矩阵中
    ferro_defense = [g for g in ferro_defense if g in adata.var_names]
    oxphos_genes = [g for g in oxphos_genes if g in adata.var_names]
    glyco_genes = [g for g in glyco_genes if g in adata.var_names]
    
    print(f"✅ 基因集过滤完毕。Ferro: {len(ferro_defense)}, OXPHOS: {len(oxphos_genes)}, Glyco: {len(glyco_genes)}")
    
    # --- 进行 Scanpy 模块打分 ---
    sc.tl.score_genes(adata, gene_list=ferro_defense, score_name='Ferro_Defense_Score')
    sc.tl.score_genes(adata, gene_list=oxphos_genes, score_name='OXPHOS_Score')
    sc.tl.score_genes(adata, gene_list=glyco_genes, score_name='Glycolysis_Score')
    
    # --- 寻找 NDUFB7 并计算机制相关性 ---
    target_gene = next((g for g in adata.var_names if g.upper() == 'NDUFB7'), None)
    
    if target_gene:
        # 提取表达量向量
        ndufb7_expr = adata[:, target_gene].X.toarray().flatten() if hasattr(adata[:, target_gene].X, 'toarray') else adata[:, target_gene].X.flatten()
        ferro_score = adata.obs['Ferro_Defense_Score'].values
        oxphos_score = adata.obs['OXPHOS_Score'].values
        glyco_score = adata.obs['Glycolysis_Score'].values
        
        # 计算 Pearson 相关系数与 p-value
        r_ferro, p_ferro = stats.pearsonr(ndufb7_expr, ferro_score)
        r_oxphos, p_oxphos = stats.pearsonr(ndufb7_expr, oxphos_score)
        r_glyco, p_glyco = stats.pearsonr(ndufb7_expr, glyco_score)
        
        print("\n==========================================")
        print("🌟 机制闭环：空间分辨率相关性结果 🌟")
        print(f" 1. NDUFB7 vs OXPHOS 产能: r = {r_oxphos:.4f} (p={p_oxphos:.2e})")
        print(f" 2. NDUFB7 vs 铁死亡防线: r = {r_ferro:.4f} (p={p_ferro:.2e})")
        print(f" 3. NDUFB7 vs 糖酵解代偿: r = {r_glyco:.4f} (p={p_glyco:.2e})")
        print("==========================================")
        
        # 逻辑判断生成结论
        if r_oxphos > 0.2 and r_ferro > 0.2:
            print("\n🎉 完美闭环！NDUFB7 空间表达与 OXPHOS 能力及铁死亡防线高度正相关。")
            print("这意味着在 NDUFB7 丢失的纤维化区，能量代谢发生了整体瘫痪，细胞的抗氧化防线也随之崩溃（易感铁死亡）。")
        else:
            print("\n⚠️ 机制结果不显著，请检查数据噪声或考虑非线性关系。")
            
    else:
        print("❌ 未找到 NDUFB7 基因。")
else:
    print("❌ 未找到 .h5ad 空间切片文件。")
