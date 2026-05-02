#!/bin/bash
# Downloads清理执行脚本
# 生成时间: $(date)
# ⚠️  执行前确认: 这些文件已成功复制到项目目录且验证通过

DOWNLOADS="$HOME/Downloads"
PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"

echo "=== Downloads清理执行 ==="
echo "预计释放空间: ~30G"
echo ""

declare -a TO_DELETE=(
    "GSE183852_DCM_Nuclei.Robj.gz"
    "GSE57338_RAW.tar"
    "GSE57338_RAW(1).tar"
    "GSE168742_RAW.tar"
    "GSE168742_human_control_CM.csv.gz"
    "GSE168742_human_HF_CM.csv.gz"
    "GSE168742_MI_day1_1.csv.gz"
    "GSE168742_MI_day1_2.csv.gz"
    "GSE168742_MI_day7_1.csv.gz"
    "GSE168742_NCM_10X_sham1.csv.gz"
    "GSE168742_NCM_10X_sham2.csv.gz"
    "GSE168742_NCM_SS2_sham1.csv.gz"
    "GSE168742_NCM_SS2_sham2.csv.gz"
    "GSE168742_NCM_SS2_TAC1.csv.gz"
    "GSE168742_NCM_SS2_TAC2.csv.gz"
    "GSE168742_WT_sham1_CM.csv.gz"
    "GSE168742_WT_sham1_FB.csv.gz"
    "GSE168742_WT_sham2_CM.csv.gz"
    "GSE168742_WT_sham2_FB.csv.gz"
    "GSE168742_WT_TAC1_CM.csv.gz"
    "GSE168742_WT_TAC1_FB.csv.gz"
    "GSE42955_RAW.tar"
    "GSE135805_RAW.tar"
    "GSE135805_series_matrix.txt.gz"
    "GSE221698_processed_data.tar.gz"
    "GSE246410_RAW(1).tar"
    "GSE246410_RAW.tar"
    "GSE246410_RAW.X-3OAAf6.tar.part"
    "GSE269054_processed_data.tar.gz"
    "pride_datasets.json"
    "722ed8b1-6bb7-4578-b887-8f264eaf01a7.h5ad"
)

released=0
for file in "${TO_DELETE[@]}"; do
    path="$DOWNLOADS/$file"
    if [ -f "$path" ]; then
        size=$(du -sb "$path" | cut -f1)
        released=$((released + size))
        rm -f "$path"
        echo "  🗑️  已删除: $file ($(du -sh "$path" 2>/dev/null | cut -f1))"
    fi
done

echo ""
echo "=== 清理完成 ==="
echo "释放空间: $(numfmt --to=iec $released 2>/dev/null || echo "$released bytes")"
echo "Downloads剩余文件:"
ls -lh "$DOWNLOADS" | tail -5
