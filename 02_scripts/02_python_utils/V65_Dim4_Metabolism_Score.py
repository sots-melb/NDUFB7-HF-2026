import scanpy as sc
import os

print("▶ 启动代谢与铁死亡联合评分...")
project_dir = os.path.expanduser("~/Projects/NDUFB7_HF_2026_04_20")
h5ad_paths = [os.path.join(dp, f) for dp, dn, filenames in os.walk(project_dir) for f in filenames if f.endswith('.h5ad') and 'Kuppe' in dp]

if h5ad_paths:
    adata = sc.read_h5ad(h5ad_paths[0])
    # 铁死亡核心5基因 (V64 确立)
    ferro_genes = ['FTL', 'SAT1', 'NFE2L2', 'SLC7A11', 'ACSL4']
    # 糖酵解与OXPHOS核心基因 (代理代谢抑制)
    glyco_genes = ['HK1', 'HK2', 'PFKM', 'ALDOA', 'ENO1', 'PKM', 'LDHA']
    oxphos_genes = ['ATP5F1A', 'COX4I1', 'CYC1', 'NDUFA4', 'SDHA']
    
    sc.tl.score_genes(adata, gene_list=[g for g in ferro_genes if g in adata.var_names], score_name='Ferro_Defense_Score')
    sc.tl.score_genes(adata, gene_list=[g for g in oxphos_genes if g in adata.var_names], score_name='OXPHOS_Score')
    
    print("✅ 代谢代理评分完成！后续可直接比较 NDUFB7-low 区域的 OXPHOS 与铁死亡防御能力是否双双崩溃。")
else:
    print("⚠️ 暂未找到 .h5ad 文件，跳过评分。")
