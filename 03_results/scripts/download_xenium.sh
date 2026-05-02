#!/bin/bash
echo "=== GSE290577 Xenium下载脚本 ==="
echo "注意: 文件可能极大(>50G)，请确认磁盘空间充足"
df -h ~/Projects/

# GEO FTP路径示例（需根据实际页面确认）
# GSE290577通常有多个supplementary files
echo ""
echo "步骤:"
echo "1. 访问 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290577"
echo "2. 找到 'Supplementary file' 区域"
echo "3. 复制下载链接（通常是ftp://ftp.ncbi.nlm.nih.gov/geo/...）"
echo "4. 用wget -c 断点续传下载"
echo ""
echo "示例命令（链接需替换为实际URL）:"
echo 'wget -c --tries=10 -P ~/Projects/NDUFB7_HF_2026_04_20/01_data/06_xenium/ "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE29nnn/GSE290577/suppl/GSE290577_xxx.tar.gz"'
