#!/bin/bash
cd ~/Projects/NDUFB7_HF_{2026_04_20}

echo "========================================"
echo "  根据筛选结果生成下载脚本"
echo "========================================"

P0_CSV="03_results/16_spatial_from_bulk/02_P0_human_heart_spatial.csv"
P1_CSV="03_results/16_spatial_from_bulk/03_P1_human_heart_sc.csv"

mkdir -p 02_scripts/10_download_jobs 01_data/04_spatial_raw 04_logs

# --- 生成P0下载脚本 ---
if [ -f "$P0_CSV" ]; then
    echo "⏳ 生成P0空间数据下载脚本..."
    
    # 提取前20个GSE
    gse_list=$(tail -n +2 "$P0_CSV" | head -20 | cut -d',' -f1 | sed 's/"//g')
    
    cat > 02_scripts/10_download_jobs/01_download_P0_spatial.sh << 'INNEREOF'
#!/bin/bash
# P0-人心脏空间转录组数据下载脚本
# 生成时间: $(date)
# 策略: 先下载Series Matrix(小文件快速验证)，再下载Supplementary

PROJECT_ROOT=~/Projects/NDUFB7_HF_{2026_04_20}
cd $PROJECT_ROOT
mkdir -p 01_data/04_spatial_raw 04_logs/spatial_download

echo "========================================"
echo "  P0空间数据批量下载"
echo "  开始时间: $(date)"
echo "========================================"

# GSE列表（从筛选结果自动生成）
GSE_LIST="
INNEREOF

    # 插入GSE列表
    echo "$gse_list" >> 02_scripts/10_download_jobs/01_download_P0_spatial.sh
    
    cat >> 02_scripts/10_download_jobs/01_download_P0_spatial.sh << 'INNEREOF'
"

# 为每个GSE创建下载函数
for gse in $GSE_LIST; do
    [ -z "$gse" ] && continue
    
    echo ""
    echo ">>> 处理 $gse <<<"
    mkdir -p "01_data/04_spatial_raw/$gse"
    
    # 1. 下载Series Matrix（小文件，含metadata）
    echo "⏳ 下载 Series Matrix..."
    nohup wget -c -t 10 --timeout=300 \
        "ftp://ftp.ncbi.nlm.nih.gov/geo/series/${gse:0:6}nnn/$gse/matrix/${gse}_series_matrix.txt.gz" \
        -O "01_data/04_spatial_raw/$gse/${gse}_series_matrix.txt.gz" \
        > "04_logs/spatial_download/${gse}_matrix.log" 2>&1 &
    
    # 2. 探测并下载Supplementary（后台，大文件）
    echo "⏳ 探测 Supplementary 文件列表..."
    # 使用GEO FTP列出suppl目录
    suppl_url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${gse:0:6}nnn/$gse/suppl/"
    
    # 获取文件列表（保存到临时文件）
    curl -s "$suppl_url" > /tmp/${gse}_ftp_list.txt 2>/dev/null
    
    if [ -s /tmp/${gse}_ftp_list.txt ]; then
        echo "   发现Supplementary文件:"
        cat /tmp/${gse}_ftp_list.txt | head -10
        
        # 优先下载小文件(<100MB)进行快速验证
        while read -r line; do
            fname=$(echo "$line" | awk '{print $NF}')
            fsize=$(echo "$line" | awk '{print $5}')
            
            # 跳过过大的文件(>2GB)先不下载，人工确认后再下
            if [ "$fsize" -gt 2147483648 ]; then
                echo "   ⏭️  $fname 过大(${fsize}bytes)，跳过待人工确认"
                continue
            fi
            
            # 优先关键词: .h5, .csv.gz, .tar.gz, .rds, .Robj
            if echo "$fname" | grep -qiE "\.(h5|csv\.gz|tar\.gz|rds|Robj|txt\.gz)$"; then
                echo "   ⬇️  下载 $fname (${fsize}bytes)..."
                nohup wget -c -t 10 --timeout=600 \
                    "${suppl_url}${fname}" \
                    -O "01_data/04_spatial_raw/$gse/${fname}" \
                    > "04_logs/spatial_download/${gse}_${fname}.log" 2>&1 &
            fi
        done < /tmp/${gse}_ftp_list.txt
    else
        echo "   ⚠️ 无法列出FTP目录，尝试HTTPS..."
        # 备选: 使用GEO页面抓取（复杂，建议人工）
    fi
    
    echo "✅ $gse 下载任务已后台启动"
    sleep 2  # 避免同时开太多连接
done

echo ""
echo "========================================"
echo "  所有P0下载任务已启动"
echo "  监控命令: tail -f 04_logs/spatial_download/*.log"
echo "  查看进度: ls -lh 01_data/04_spatial_raw/*/"
echo "========================================"
INNEREOF

    chmod +x 02_scripts/10_download_jobs/01_download_P0_spatial.sh
    echo "✅ P0下载脚本已生成: 02_scripts/10_download_jobs/01_download_P0_spatial.sh"
    echo "   包含GSE: $(echo "$gse_list" | wc -l) 个"
else
    echo "❌ P0 CSV不存在，跳过生成"
fi

# --- 生成P1下载脚本（简化版） ---
if [ -f "$P1_CSV" ]; then
    echo "⏳ 生成P1单细胞下载脚本..."
    gse_list_p1=$(tail -n +2 "$P1_CSV" | head -10 | cut -d',' -f1 | sed 's/"//g')
    
    cat > 02_scripts/10_download_jobs/02_download_P1_sc.sh << 'INNEREOF'
#!/bin/bash
# P1-人心脏单细胞数据下载脚本
cd ~/Projects/NDUFB7_HF_{2026_04_20}
mkdir -p 01_data/05_sc_reference 04_logs/sc_download

GSE_LIST="
INNEREOF
    echo "$gse_list_p1" >> 02_scripts/10_download_jobs/02_download_P1_sc.sh
    cat >> 02_scripts/10_download_jobs/02_download_P1_sc.sh << 'INNEREOF'
"

for gse in $GSE_LIST; do
    [ -z "$gse" ] && continue
    echo "⏳ 下载 $gse Series Matrix..."
    mkdir -p "01_data/05_sc_reference/$gse"
    nohup wget -c -t 10 --timeout=300 \
        "ftp://ftp.ncbi.nlm.nih.gov/geo/series/${gse:0:6}nnn/$gse/matrix/${gse}_series_matrix.txt.gz" \
        -O "01_data/05_sc_reference/$gse/${gse}_series_matrix.txt.gz" \
        > "04_logs/sc_download/${gse}.log" 2>&1 &
    sleep 1
done
echo "✅ P1下载已启动"
INNEREOF

    chmod +x 02_scripts/10_download_jobs/02_download_P1_sc.sh
    echo "✅ P1下载脚本已生成"
fi

echo ""
echo "========================================"
echo "  下载脚本生成完成"
echo "========================================"
echo ""
echo "【使用指南】"
echo "  1. 先人工确认P0 Top 3的GEO页面（浏览器检查Supplementary）"
echo "  2. 确认无误后执行: bash 02_scripts/10_download_jobs/01_download_P0_spatial.sh"
echo "  3. 所有下载后台运行，不影响前台工作"
echo "  4. 30分钟后运行QC脚本验证"
echo ""
echo "【监控命令】"
echo "  tail -f 04_logs/spatial_download/*.log"
echo "  watch -n 30 'ls -lh 01_data/04_spatial_raw/*/'"
