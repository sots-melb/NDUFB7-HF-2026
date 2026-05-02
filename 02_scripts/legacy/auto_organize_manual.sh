#!/bin/bash
PROJECT_ROOT="$HOME/Projects/NDUFB7_HF_{2026_04_20}"
DOWNLOAD_DIR="$HOME/Downloads"
DATA_RAW="$PROJECT_ROOT/01_data/01_raw_geo"

echo "=== 手动下载文件归集 ==="
echo "扫描目录: $DOWNLOAD_DIR"

# 移动GSE141910相关文件
if ls $DOWNLOAD_DIR/*141910* 1> /dev/null 2>&1; then
    mkdir -p "$DATA_RAW/GSE141910"
    mv -v $DOWNLOAD_DIR/*141910* "$DATA_RAW/GSE141910/"
    echo "[归集] GSE141910 文件已移动"
fi

# 移动GSE57338相关文件
if ls $DOWNLOAD_DIR/*57338* 1> /dev/null 2>&1; then
    mkdir -p "$DATA_RAW/GSE57338"
    mv -v $DOWNLOAD_DIR/*57338* "$DATA_RAW/GSE57338/"
    echo "[归集] GSE57338 文件已移动"
fi

# 移动GSE5406相关文件
if ls $DOWNLOAD_DIR/*5406* 1> /dev/null 2>&1; then
    mkdir -p "$DATA_RAW/GSE5406"
    mv -v $DOWNLOAD_DIR/*5406* "$DATA_RAW/GSE5406/"
    echo "[归集] GSE5406 文件已移动"
fi

# 显示当前资产状态
echo ""
echo "=== 当前数据资产状态 ==="
for d in GSE141910 GSE57338 GSE116250 GSE5406 GSE46224 GSE48166 GSE55296 GSE79962; do
    dir_path="$DATA_RAW/$d"
    if [ -d "$dir_path" ]; then
        size=$(du -sh "$dir_path" 2>/dev/null | cut -f1)
        file_count=$(ls -1 "$dir_path" 2>/dev/null | wc -l)
        echo "[$d] $size ($file_count 个文件)"
        ls -lh "$dir_path" | tail -n +2 | awk '{print "  - " $9 " (" $5 ")"}'
    fi
done
