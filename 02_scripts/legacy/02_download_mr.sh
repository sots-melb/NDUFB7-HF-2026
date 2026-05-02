
#!/bin/bash
# MR数据下载脚本（手动执行或放入后台）
# 建议逐个下载，避免网络拥堵

mkdir -p ~/Projects/NDUFB7_HF_{2026_04_20}/01_data/03_mr_data

# 1. eQTLGen (需从网站手动下载，无直接FTP)
# 访问: https://www.eqtlgen.org/phase1.html
# 下载: cis-eQTL summary statistics
# 搜索基因: NDUFB7 (chr19: 需确认具体位置)

# 2. GTEx心脏eQTL (可直接wget)
cd ~/Projects/NDUFB7_HF_{2026_04_20}/01_data/03_mr_data
wget -c https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_all_associations/Heart_Left_Ventricle.v8.allpairs.txt.gz
echo "GTEx心脏eQTL下载完成"

# 3. FinnGen心衰GWAS (需从网站下载，或申请API)
# 访问: https://www.finngen.fi/en/access_results
# 下载: Endpoint: I9_HEARTFAIL (Heart failure)

# 4. 备选: UK Biobank GWAS summary
# 需通过GWAS Catalog或Neale Lab获取

echo "MR数据下载脚本准备完成"

