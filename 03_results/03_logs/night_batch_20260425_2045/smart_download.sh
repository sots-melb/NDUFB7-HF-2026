#!/bin/bash
# 智能断点续传下载函数
# 用法: smart_download <URL> <输出文件> <日志文件> [最大重试次数]
smart_download() {
    local url="$1"
    local output="$2"
    local log="$3"
    local max_try="${4:-20}"
    local count=0
    local wait_time=120  # 失败后等待2分钟
    
    echo "[$(date '+%H:%M:%S')] START: $output from $url" >> "$log"
    
    while [ $count -lt $max_try ]; do
        local attempt=$((count+1))
        echo "[$(date '+%H:%M:%S')] Attempt $attempt/$max_try: downloading..." >> "$log"
        
        # wget断点续传：-c继续，--timeout=60防止挂死，--tries=1配合循环重试
        wget --timeout=120 --tries=1 --continue "$url" -O "$output" 2>> "$log"
        local exit_code=$?
        
        # 检查是否成功（文件存在且非空）
        if [ $exit_code -eq 0 ] && [ -s "$output" ]; then
            local sz=$(du -sh "$output" 2>/dev/null | cut -f1)
            echo "[$(date '+%H:%M:%S')] SUCCESS: $output ($sz)" >> "$log"
            return 0
        else
            echo "[$(date '+%H:%M:%S')] FAILED (exit=$exit_code, size=$(du -sh $output 2>/dev/null | cut -f1)), retry in ${wait_time}s..." >> "$log"
            sleep $wait_time
            count=$((count+1))
            # 渐进式增加等待时间（网络持续故障时）
            if [ $count -gt 5 ]; then wait_time=300; fi
            if [ $count -gt 10 ]; then wait_time=600; fi
        fi
    done
    
    echo "[$(date '+%H:%M:%S')] GAVE UP: $output after $max_try attempts" >> "$log"
    return 1
}
