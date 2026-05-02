#!/bin/bash
PROJECT="/home/y411869/Projects/NDUFB7_HF_2026_04_20"
echo "========================================"
echo "V171C: GSE157282 资产定位"
echo "========================================"

# 全面搜索
echo "--- 在Project中搜索 ---"
find "$PROJECT" -name "*157282*" 2>/dev/null

echo "--- 在Downloads中搜索 ---"
find ~/Downloads -name "*157282*" 2>/dev/null

echo "--- 检查GEO标准路径 ---"
ls -la "$PROJECT/01_data/01_raw_geo/GSE157282/" 2>/dev/null || echo "  Directory not found"

echo "--- 检查soft文件 ---"
find "$PROJECT" ~/Downloads -name "*157282*soft*" -o -name "*157282*matrix*" 2>/dev/null

echo ""
echo "如果以上都没有，需要重新下载："
echo "  手动: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157282"
echo "  或运行: cd ~/Downloads && wget -c https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE157282&format=file -O GSE157282_RAW.tar"
