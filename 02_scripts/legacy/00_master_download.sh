#!/bin/bash
# NDUFB7项目数据资产总控脚本 v1.0
# 功能：批量下载GEO、断点续传、完整性验证、自动归档、R脚本生成、环境快照
# 助教提示：此脚本只需执行一次，中断后可重复执行，不会重复下载已完成的文件

set -euo pipefail

# ==================== 配置区 ====================
PROJECT_ROOT="$HOME/Projects/NDUFB7_HF_{2026_04_20}"
LOG_DIR="$PROJECT_ROOT/05_logs"
DATA_RAW="$PROJECT_ROOT/01_data/01_raw_geo"
SCRIPTS_QC="$PROJECT_ROOT/02_scripts/02_qc"
SCRIPTS_ANALYSIS="$PROJECT_ROOT/02_scripts/03_analysis"
TMP_DIR="$PROJECT_ROOT/06_tmp"

mkdir -p "$LOG_DIR" "$DATA_RAW" "$SCRIPTS_QC" "$SCRIPTS_ANALYSIS" "$TMP_DIR"

# 日志文件：每次执行生成独立日志，永不覆盖
DOWNLOAD_LOG="$LOG_DIR/batch_download_$(date +%Y%m%d_%H%M%S).log"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] === NDUFB7数据资产总控脚本启动 ===" | tee -a "$DOWNLOAD_LOG"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] 项目根目录: $PROJECT_ROOT" | tee -a "$DOWNLOAD_LOG"

# 下载工具检测（优先wget，备选curl）
DL_TOOL=""
if command -v wget &> /dev/null; then
    DL_TOOL="wget"
    echo "[环境] 使用 wget 下载" | tee -a "$DOWNLOAD_LOG"
elif command -v curl &> /dev/null; then
    DL_TOOL="curl"
    echo "[环境] 使用 curl 下载" | tee -a "$DOWNLOAD_LOG"
else
    echo "[致命错误] 系统未安装wget或curl，请联系服务器管理员" | tee -a "$DOWNLOAD_LOG"
    exit 1
fi

# GEO FTP路径计算函数：自动将GSE编号转换为正确的FTP子目录
# 规则：GSE数字补零到6位，取前3位+nnn，例如GSE57338 -> 057338 -> GSE057nnn
get_geo_ftp_dir() {
    local gse=$1
    local num=${gse#GSE}
    local num_len=${#num}
  local prefix=${num:0:3}
    local prefix=${num:0:3}
    echo "GSE${prefix}nnn"
}

# 下载函数：Series Matrix（断点续传+完整性验证）
download_series_matrix() {
    local gse=$1
    local ftp_dir=$(get_geo_ftp_dir "$gse")
    local ftp_url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${ftp_dir}/${gse}/matrix/${gse}_series_matrix.txt.gz"
    local target_dir="$DATA_RAW/$gse"
    local target_file="$target_dir/${gse}_series_matrix.txt.gz"
    
    mkdir -p "$target_dir"
    
    # 检查是否已完整下载（通过gunzip -t验证）
    if [ -f "$target_file" ]; then
        if gunzip -t "$target_file" 2>/dev/null; then
            local local_size=$(stat -c%s "$target_file" 2>/dev/null || stat -f%z "$target_file" 2>/dev/null)
            echo "[跳过] $gse 已存在且校验通过 (${local_size} bytes)" | tee -a "$DOWNLOAD_LOG"
            return 0
        else
            echo "[断点续传] $gse 文件损坏/不完整，重新下载..." | tee -a "$DOWNLOAD_LOG"
            rm -f "$target_file"
        fi
    fi
    
    echo "[下载] $gse Series Matrix..." | tee -a "$DOWNLOAD_LOG"
    echo "  URL: $ftp_url" | tee -a "$DOWNLOAD_LOG"
    
    # 执行下载：-c断点续传, -t 0无限重试, --timeout=60防假死
    if [ "$DL_TOOL" == "wget" ]; then
        wget -c -t 0 --timeout=60 --progress=dot:giga -O "$target_file.tmp" "$ftp_url" 2>>"$DOWNLOAD_LOG"
    else
        curl -C - --retry 0 --max-time 60 -o "$target_file.tmp" "$ftp_url" 2>>"$DOWNLOAD_LOG"
    fi
    
    if [ $? -eq 0 ] && [ -f "$target_file.tmp" ]; then
        mv "$target_file.tmp" "$target_file"
        local final_size=$(stat -c%s "$target_file" 2>/dev/null || stat -f%z "$target_file" 2>/dev/null)
        echo "[完成] $gse 下载成功 (${final_size} bytes)" | tee -a "$DOWNLOAD_LOG"
        
        # 完整性验证：gunzip -t 测试gzip文件完整性
        if gunzip -t "$target_file" 2>/dev/null; then
            echo "[验证] $gse gzip完整性校验通过" | tee -a "$DOWNLOAD_LOG"
            # 记录MD5到注册表
            md5sum "$target_file" >> "$LOG_DIR/md5_registry.txt"
            return 0
        else
            echo "[错误] $gse 文件校验失败，已删除，请重新执行脚本" | tee -a "$DOWNLOAD_LOG"
            rm -f "$target_file"
            return 1
        fi
    else
        echo "[错误] $gse 下载失败，网络或GEO服务器问题" | tee -a "$DOWNLOAD_LOG"
        rm -f "$target_file.tmp"
        return 1
    fi
}

# 单细胞/特殊数据下载：Supplementary processed files
download_supplementary() {
    local gse=$1
    local ftp_dir=$(get_geo_ftp_dir "$gse")
    local suppl_url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${ftp_dir}/${gse}/suppl/"
    local target_dir="$DATA_RAW/$gse"
    
    mkdir -p "$target_dir"
    echo "[下载] $gse 检查Supplementary files..." | tee -a "$DOWNLOAD_LOG"
    
    # 获取文件列表
    local file_list="$TMP_DIR/${gse}_suppl_list.txt"
    if [ "$DL_TOOL" == "wget" ]; then
        wget -q --timeout=60 -O - "$suppl_url" 2>/dev/null | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"$//' > "$file_list"
    else
        curl -s --max-time 60 "$suppl_url" 2>/dev/null | grep -oE 'href="[^"]+"' | sed 's/href="//;s/"$//' > "$file_list"
    fi
    
    if [ ! -s "$file_list" ]; then
        echo "[警告] $gse 无法获取补充文件列表（可能无suppl或网络延迟）" | tee -a "$DOWNLOAD_LOG"
        return 1
    fi
    
    echo "[信息] $gse 发现补充文件：" | tee -a "$DOWNLOAD_LOG"
    cat "$file_list" | tee -a "$DOWNLOAD_LOG"
    
    # 智能匹配：根据GSE特性下载关键processed文件
    local files_to_download=()
    case "$gse" in
        GSE168742)
            # 单细胞processed CSV（人类raw不可用，但processed可用）
            files_to_download=("human_HF_CM.csv.gz" "human_control_CM.csv.gz" "mouse_CM.csv.gz" "mouse_TAC_CM.csv.gz")
            ;;
        GSE315590)
            # 时间序列metadata
            files_to_download=("metadata.txt.gz" "TAC_transition")
            ;;
        GSE165303)
            # RNA-seq master table
            files_to_download=("RNA-Seq_master_table.txt.gz" "counts.txt.gz")
            ;;
        GSE270788)
            files_to_download=("Multiome" "RNA_ATAC" "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz")
            ;;
    esac
    
    for pattern in "${files_to_download[@]}"; do
        local matched_file=$(grep "$pattern" "$file_list" | head -1)
        if [ -n "$matched_file" ]; then
            local f_url="${suppl_url}${matched_file}"
            local f_path="$target_dir/$matched_file"
            
            if [ -f "$f_path" ]; then
                echo "[跳过] $matched_file 已存在" | tee -a "$DOWNLOAD_LOG"
            else
                echo "[下载] 补充文件: $matched_file ..." | tee -a "$DOWNLOAD_LOG"
                if [ "$DL_TOOL" == "wget" ]; then
                    wget -c -t 0 --timeout=60 --progress=dot:giga -O "$f_path.tmp" "$f_url" 2>>"$DOWNLOAD_LOG" && mv "$f_path.tmp" "$f_path"
                else
                    curl -C - --max-time 60 -o "$f_path.tmp" "$f_url" 2>>"$DOWNLOAD_LOG" && mv "$f_path.tmp" "$f_path"
                fi
                
                if [ -f "$f_path" ]; then
                    echo "[完成] $matched_file 下载成功" | tee -a "$DOWNLOAD_LOG"
                    md5sum "$f_path" >> "$LOG_DIR/md5_registry.txt" 2>/dev/null || true
                fi
            fi
        fi
    done
}

# 迁移旧数据：自动将~/Downloads中的相关文件归档到项目区
migrate_downloads() {
    local download_dir="$HOME/Downloads"
    if [ ! -d "$download_dir" ]; then
        return 0
    fi
    
    echo "[迁移] 扫描 ~/Downloads 中的GEO文件..." | tee -a "$DOWNLOAD_LOG"
    local all_gse=("${BULK_GSE_LIST[@]}" "${SC_GSE_LIST[@]}")
    for gse in "${all_gse[@]}"; do
        local found=$(find "$download_dir" -maxdepth 1 -type f -iname "*${gse}*" 2>/dev/null)
        if [ -n "$found" ]; then
            mkdir -p "$DATA_RAW/$gse"
            echo "$found" | while read -r f; do
                local fname=$(basename "$f")
                if [ ! -f "$DATA_RAW/$gse/$fname" ]; then
                    mv "$f" "$DATA_RAW/$gse/"
                    echo "[迁移] $fname -> $DATA_RAW/$gse/" | tee -a "$DOWNLOAD_LOG"
                fi
            done
        fi
    done
}

# ==================== 主执行区 ====================

# 1. 先迁移旧数据（如果你之前用火狐下载过）
migrate_downloads

# 2. 定义下载队列
BULK_GSE_LIST=("GSE141910" "GSE57338" "GSE116250" "GSE5406" "GSE46224" "GSE48166" "GSE55296" "GSE79962")
SC_GSE_LIST=("GSE168742" "GSE315590" "GSE165303" "GSE270788")

# 3. 批量下载Bulk数据集
echo "" | tee -a "$DOWNLOAD_LOG"
echo "========== Phase 1: 下载核心Bulk数据集 ==========" | tee -a "$DOWNLOAD_LOG"
for gse in "${BULK_GSE_LIST[@]}"; do
    download_series_matrix "$gse"
    sleep 2  # 礼貌间隔，避免NCBI服务器封IP
done

# 4. 批量下载单细胞/特殊数据集
echo "" | tee -a "$DOWNLOAD_LOG"
echo "========== Phase 2: 下载单细胞/特殊数据集 ==========" | tee -a "$DOWNLOAD_LOG"
for gse in "${SC_GSE_LIST[@]}"; do
    # 单细胞可能没有Series Matrix，不中断流程
    download_series_matrix "$gse" || true
    download_supplementary "$gse"
    sleep 2
done

# ==================== 数据资产报告 ====================
echo "" | tee -a "$DOWNLOAD_LOG"
echo "========== Phase 3: 生成数据资产报告 ==========" | tee -a "$DOWNLOAD_LOG"

REPORT_FILE="$LOG_DIR/data_asset_report_$(date +%Y%m%d_%H%M%S).txt"
{
    echo "NDUFB7项目数据资产报告"
    echo "生成时间: $(date)"
    echo "操作者: $(whoami)"
    echo "====================================="
    echo ""
    echo "【下载日志】"
    cat "$DOWNLOAD_LOG"
    echo ""
    echo "【文件清单与大小】"
    for gse_dir in "$DATA_RAW"/GSE*; do
        if [ -d "$gse_dir" ]; then
            echo ""
            echo "[$(basename "$gse_dir")]"
            ls -lh "$gse_dir" 2>/dev/null || echo "  (空目录)"
        fi
    done
    echo ""
    echo "【MD5注册表】"
    cat "$LOG_DIR/md5_registry.txt" 2>/dev/null || echo "暂无"
} > "$REPORT_FILE"

echo "[报告] 完整资产报告已生成: $REPORT_FILE" | tee -a "$DOWNLOAD_LOG"

# ==================== R脚本自动生成 ====================
echo "" | tee -a "$DOWNLOAD_LOG"
echo "========== Phase 4: 生成标准化R脚本 ==========" | tee -a "$DOWNLOAD_LOG"

# R脚本1: Bulk数据导入与QC（保存到02_scripts/02_qc/）
cat > "$SCRIPTS_QC/01_bulk_data_import_qc_v1.R" << 'REOF'
# ==============================================================================
# 项目: NDUFB7心衰多组学研究
# 脚本: 01_bulk_data_import_qc_v1.R
# 阶段: Phase 1 - 数据导入与QC
# 功能: 批量读取GEO Series Matrix，提取表达矩阵和表型，基础QC与NDUFB7检测
# 生成日期: $(date +%Y-%m-%d)
# 状态: 已就绪
# 上游: 00_master_download.sh (数据下载)
# 下游: 02_preprocess_bulk_v1.R (标准化)
# ==============================================================================

PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

# 结果目录
RESULT_DIR <- file.path(PROJECT_ROOT, "03_results/01_qc_reports")
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

# 日志：同时输出到控制台和文件
LOG_FILE <- file.path(PROJECT_ROOT, "05_logs", paste0("R_bulkQC_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(LOG_FILE, split = TRUE)
cat("[", format(Sys.time()), "] 脚本启动: 01_bulk_data_import_qc_v1.R\n")
cat("[", format(Sys.time()), "] 工作目录:", getwd(), "\n")

# 加载包（失败时给出友好提示）
cat("[", format(Sys.time()), "] 加载R包...\n")
pkg_list <- c("GEOquery", "limma", "tidyverse", "ggplot2")
for (pkg in pkg_list) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("[错误] 缺少包:", pkg, "请执行 BiocManager::install('", pkg, "')\n")
  } else {
    library(pkg, character.only = TRUE)
    cat("[OK] 已加载:", pkg, "\n")
  }
}

RAW_DIR <- file.path(PROJECT_ROOT, "01_data/01_raw_geo")
GSE_BULK <- c("GSE141910", "GSE57338", "GSE116250", "GSE5406", "GSE46224", "GSE48166", "GSE55296", "GSE79962")

geo_list <- list()

for (gse_id in GSE_BULK) {
  cat("\n========== [处理]", gse_id, "==========\n")
  sm_file <- file.path(RAW_DIR, gse_id, paste0(gse_id, "_series_matrix.txt.gz"))
  
  if (!file.exists(sm_file)) {
    cat("[跳过] 文件不存在:", sm_file, "\n")
    next
  }
  
  # 读取GEO Series Matrix
  gse <- tryCatch({
    getGEO(filename = sm_file, GSEMatrix = TRUE)
  }, error = function(e) {
    cat("[错误] getGEO失败:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(gse) || length(gse) == 0) {
    cat("[跳过] 读取结果为空\n")
    next
  }
  
  exprs <- exprs(gse[[1]])
  pheno <- pData(gse[[1]])
  
  cat("[信息] 表达矩阵:", nrow(exprs), "基因 x", ncol(exprs), "样本\n")
  cat("[信息] 表型字段:", paste(head(colnames(pheno), 5), collapse = ", "), "...\n")
  
  # 保存解析后的RDS（下游脚本直接读取此文件，无需重复解析GEO）
  rds_path <- file.path(RAW_DIR, gse_id, paste0(gse_id, "_parsed.rds"))
  saveRDS(list(exprs = exprs, pheno = pheno, gse = gse[[1]]), file = rds_path)
  cat("[保存] 已导出:", rds_path, "\n")
  
  # NDUFB7检测（核心QC）
  ndufb7_probe <- grep("^NDUFB7$", rownames(exprs), value = TRUE, ignore.case = TRUE)
  if (length(ndufb7_probe) > 0) {
    ndufb7_expr <- as.numeric(exprs[ndufb7_probe[1], ])
    cat("[QC-通过] NDUFB7检测成功! 表达量范围:", round(range(ndufb7_expr, na.rm = TRUE), 3), "\n")
    
    # 快速箱线图保存到QC结果目录
    pdf(file.path(RESULT_DIR, paste0(gse_id, "_NDUFB7_boxplot.pdf")), width = 6, height = 4)
    boxplot(ndufb7_expr ~ pheno$title, main = paste(gse_id, "NDUFB7 Expression"), las = 2, cex.axis = 0.7)
    dev.off()
  } else {
    cat("[QC-警告] NDUFB7未在行名中找到，可能需要探针ID转换（尤其GSE5406的U133A芯片）\n")
  }
  
  geo_list[[gse_id]] <- gse[[1]]
}

cat("\n========== [汇总] ==========\n")
cat("成功导入数据集:", length(geo_list), "/", length(GSE_BULK), "\n")
cat("已保存RDS位置:", RAW_DIR, "\n")
cat("QC图表位置:", RESULT_DIR, "\n")
cat("[下一步] 请运行 02_scripts/03_analysis/01_meta_analysis_setup_v1.R 开始荟萃分析准备\n")

sink()
cat("\n[完成] R脚本执行结束，完整日志见:", LOG_FILE, "\n")
REOF

echo "[R脚本] 已生成: $SCRIPTS_QC/01_bulk_data_import_qc_v1.R" | tee -a "$DOWNLOAD_LOG"

# R脚本2: 单细胞数据读取模板（保存到02_scripts/03_analysis/）
cat > "$SCRIPTS_ANALYSIS/01_scRNA_import_template_v1.R" << 'REOF'
# ==============================================================================
# 脚本: 01_scRNA_import_template_v1.R
# 阶段: Phase 2 - 单细胞分析准备
# 功能: GSE168742 processed CSV读取与Seurat对象创建
# 注意: GSE168742人类raw data不可用，使用processed CSV
# ==============================================================================

PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

library(Seurat)
library(tidyverse)

SC_DIR <- file.path(PROJECT_ROOT, "01_data/01_raw_geo/GSE168742")

# 读取processed CSV（如果下载成功）
hf_file <- file.path(SC_DIR, "GSE168742_human_HF_CM.csv.gz")
ctrl_file <- file.path(SC_DIR, "GSE168742_human_control_CM.csv.gz")

if (file.exists(hf_file) && file.exists(ctrl_file)) {
  cat("[信息] 读取单细胞processed CSV...\n")
  hf_data <- read_csv(hf_file)
  ctrl_data <- read_csv(ctrl_file)
  cat("[完成] HF样本:", ncol(hf_data), "细胞; Control样本:", ncol(ctrl_data), "细胞\n")
  cat("[提示] 下一步: CreateSeuratObject -> QC -> Harmony整合\n")
} else {
  cat("[警告] CSV文件未找到，请检查GSE168742 supplementary下载是否完成\n")
  cat("[提示] 如GEO无自动下载，请手动从GEO页面下载human_HF_CM.csv.gz并放入:", SC_DIR, "\n")
}
REOF

echo "[R脚本] 已生成: $SCRIPTS_ANALYSIS/01_scRNA_import_template_v1.R" | tee -a "$DOWNLOAD_LOG"

# ==================== 环境恢复脚本生成 ====================
echo "" | tee -a "$DOWNLOAD_LOG"
echo "========== Phase 5: 生成环境恢复工具 ==========" | tee -a "$DOWNLOAD_LOG"

# Bash环境恢复
cat > "$PROJECT_ROOT/load_project_env.sh" << 'EOF'
#!/bin/bash
# NDUFB7项目环境恢复脚本
# 用法: 在新Terminal中执行: source ~/Projects/NDUFB7_HF_{2026_04_20}/load_project_env.sh

PROJECT_ROOT="$HOME/Projects/NDUFB7_HF_{2026_04_20}"
export PROJECT_ROOT

cd "$PROJECT_ROOT" || exit 1
echo ""
echo "╔══════════════════════════════════════════╗"
echo "║     NDUFB7项目环境恢复成功               ║"
echo "╚══════════════════════════════════════════╝"
echo "当前目录: $(pwd)"
echo ""

# 显示数据资产状态
echo "【数据资产状态】"
if [ -d "$PROJECT_ROOT/01_data/01_raw_geo" ]; then
    du -sh "$PROJECT_ROOT/01_data/01_raw_geo"/* 2>/dev/null | sed 's/^/  /'
else
    echo "  (暂无数据)"
fi

echo ""
echo "【最近操作日志】"
ls -lt "$PROJECT_ROOT/05_logs" 2>/dev/null | head -3 | sed 's/^/  /'

echo ""
echo "【可执行的R脚本】"
find "$PROJECT_ROOT/02_scripts" -name "*.R" -type f 2>/dev/null | sed 's|.*/||' | sed 's/^/  /'

echo ""
echo "【助教提示】"
echo "  1. 在RStudio中恢复项目: File -> Open Project -> 选择 NDUFB7_HF.Rproj"
echo "  2. 运行Bulk QC: Rscript $PROJECT_ROOT/02_scripts/02_qc/01_bulk_data_import_qc_v1.R"
echo "  3. 查看日志: cat $PROJECT_ROOT/05_logs/最新日志文件"
echo ""
EOF

chmod +x "$PROJECT_ROOT/load_project_env.sh"
echo "[环境] 已生成Bash恢复脚本: $PROJECT_ROOT/load_project_env.sh" | tee -a "$DOWNLOAD_LOG"

# RStudio Project文件（双击即可恢复完整R环境）
cat > "$PROJECT_ROOT/NDUFB7_HF.Rproj" << 'EOF'
Version: 1.0

RestoreWorkspace: No
SaveWorkspace: No
AlwaysSaveHistory: Default

EnableCodeIndexing: Yes
UseSpacesForTab: Yes
NumSpacesForTab: 2
Encoding: UTF-8

RnwWeave: Sweave
LaTeX: pdfLaTeX

AutoAppendNewline: Yes
StripTrailingWhitespace: Yes
EOF

echo "[环境] 已生成RStudio项目文件: $PROJECT_ROOT/NDUFB7_HF.Rproj" | tee -a "$DOWNLOAD_LOG"

# ==================== 最终提示 ====================
echo "" | tee -a "$DOWNLOAD_LOG"
echo "╔════════════════════════════════════════════════════════════╗" | tee -a "$DOWNLOAD_LOG"
echo "║              总控脚本执行完毕                              ║" | tee -a "$DOWNLOAD_LOG"
echo "╠════════════════════════════════════════════════════════════╣" | tee -a "$DOWNLOAD_LOG"
echo "║  数据位置: $DATA_RAW" | tee -a "$DOWNLOAD_LOG"
echo "║  日志位置: $LOG_DIR" | tee -a "$DOWNLOAD_LOG"
echo "║  R脚本位置: $SCRIPTS_QC 和 $SCRIPTS_ANALYSIS" | tee -a "$DOWNLOAD_LOG"
echo "║  恢复脚本: source $PROJECT_ROOT/load_project_env.sh" | tee -a "$DOWNLOAD_LOG"
echo "╚════════════════════════════════════════════════════════════╝" | tee -a "$DOWNLOAD_LOG"

