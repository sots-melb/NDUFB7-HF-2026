#!/bin/bash
cd ~/Projects/NDUFB7_HF_{2026_04_20}

echo "╔════════════════════════════════════════════════════════════╗"
echo "║     NDUFB7空间转录组筛选流水线 - 一键执行                ║"
echo "╚════════════════════════════════════════════════════════════╝"

echo ""
echo "【阶段1/4】解析2万+GEO摘要并筛选空间候选"
echo "========================================"
bash -c 'Rscript 02_scripts/09_geo_bulk_screen/01_parse_and_screen_20k.R 2>&1 | tee 04_logs/30_parse_screen_20k.log'

echo ""
echo "【阶段2/4】生成下载脚本"
echo "========================================"
bash 02_scripts/09_geo_bulk_screen/02_generate_download_scripts.sh 2>&1 | tee 04_logs/31_generate_download.log

echo ""
echo "========================================"
echo "  流水线阶段1-2完成"
echo "========================================"
echo ""
echo "【人工检查点】"
echo "  1. 查看P0候选: cat 03_results/16_spatial_from_bulk/07_screening_report.txt"
echo "  2. 浏览器验证Top 3 GEO页面的Supplementary文件"
echo "  3. 确认无误后执行: bash 02_scripts/10_download_jobs/01_download_P0_spatial.sh"
echo ""
echo "【阶段3/4】下载完成后执行QC（自动或手动）"
echo "  Rscript 02_scripts/09_geo_bulk_screen/03_qc_verify.R"
echo ""
echo "【阶段4/4】读取可用数据并可视化NDUFB7"
echo "  Rscript 02_scripts/09_geo_bulk_screen/04_auto_read_spatial.R"
