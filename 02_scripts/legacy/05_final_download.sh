#!/bin/bash
# NDUFB7_Mito_2026 - 最终版批量下载脚本
set -e
PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"
GEO="$PROJECT/01_data/01_raw_geo"
MR="$PROJECT/01_data/03_mr_data"
LOG="$PROJECT/04_logs"
mkdir -p "$GEO" "$MR" "$LOG"

G='\033[32m'; R='\033[31m'; Y='\033[33m'; C='\033[36m'; NC='\033[0m'
info() { echo -e "${C}[$(date +%H:%M:%S)]${NC} $*"; }
ok()   { echo -e "${G}[$(date +%H:%M:%S)] ✅${NC} $*"; }
warn() { echo -e "${Y}[$(date +%H:%M:%S)] ⚠️${NC} $*"; }
err()  { echo -e "${R}[$(date +%H:%M:%S)] ❌${NC} $*"; }

download_series_matrix() {
    local gse=$1
    local dir="$GEO/$gse"
    mkdir -p "$dir"
    if [ -f "$dir/${gse}_series_matrix.txt.gz" ] && [ -s "$dir/${gse}_series_matrix.txt.gz" ]; then
        local sz=$(du -sh "$dir/${gse}_series_matrix.txt.gz" | cut -f1)
        ok "$gse 已存在 ($sz) → 跳过"
        return 0
    fi
    local num=${gse#GSE}
    local prefix="${num:0:-3}nnn"
    local url="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${gse}/matrix/${gse}_series_matrix.txt.gz"
    info "⬇️  下载 $gse ..."
    if wget -q --timeout=90 --tries=2 -c -O "$dir/${gse}_series_matrix.txt.gz" "$url" 2>>"$LOG/${gse}.log"; then
        local sz=$(du -sh "$dir/${gse}_series_matrix.txt.gz" | cut -f1)
        ok "$gse 完成 ($sz)"
        return 0
    else
        err "$gse 失败"
        return 1
    fi
}

download_soft() {
    local gse=$1
    local dir="$GEO/$gse"
    [ -f "$dir/${gse}_family.soft.gz" ] && return 0
    local num=${gse#GSE}
    local prefix="${num:0:-3}nnn"
    wget -q --timeout=60 -c -O "$dir/${gse}_family.soft.gz" \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${gse}/soft/${gse}_family.soft.gz" 2>>"$LOG/${gse}.log" || true
}

echo ""
echo -e "${C}╔══════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${C}║     NDUFB7_Mito_2026 - 48小时高速GEO批量下载                     ║${NC}"
echo -e "${C}║     开始: $(date '+%Y-%m-%d %H:%M:%S')                                    ║${NC}"
echo -e "${C}║     磁盘: $(df -h $PROJECT | awk 'NR==2{print $4}') 可用                                  ║${NC}"
echo -e "${C}╚══════════════════════════════════════════════════════════════════╝${NC}"
echo ""
START=$(date +%s)
SUCCESS=0; FAILED=0

# Phase 1: P0核心
info "═══════════════════════════════════════════════════════════════"
info "  Phase 1/5: P0核心数据集 (最高优先级)"
info "═══════════════════════════════════════════════════════════════"
for gse in GSE168742 GSE141910 GSE57338 GSE116250 GSE5406 GSE79962; do
    if download_series_matrix "$gse"; then SUCCESS=$((SUCCESS+1)); else FAILED=$((FAILED+1)); fi
    download_soft "$gse"
    sleep 2
done

# Phase 2: P1空间
info "═══════════════════════════════════════════════════════════════"
info "  Phase 2/5: P1空间转录组"
info "═══════════════════════════════════════════════════════════════"
for gse in GSE225295 GSE210867 GSE277721 GSE190132 GSE198699 GSE184494; do
    AVAIL=$(df -m "$PROJECT" | awk 'NR==2{print $4}')
    [ "$AVAIL" -lt 2048 ] && { warn "磁盘不足，暂停"; break; }
    if download_series_matrix "$gse"; then SUCCESS=$((SUCCESS+1)); else FAILED=$((FAILED+1)); fi
    download_soft "$gse"
    sleep 2
done

# Phase 3: P1单细胞
info "═══════════════════════════════════════════════════════════════"
info "  Phase 3/5: P1单细胞扩展"
info "═══════════════════════════════════════════════════════════════"
for gse in GSE183852 GSE120064 GSE165838 GSE315590 GSE46224 GSE48166 GSE217494 GSE270788; do
    AVAIL=$(df -m "$PROJECT" | awk 'NR==2{print $4}')
    [ "$AVAIL" -lt 2048 ] && { warn "磁盘不足"; break; }
    if download_series_matrix "$gse"; then SUCCESS=$((SUCCESS+1)); else FAILED=$((FAILED+1)); fi
    download_soft "$gse"
    sleep 2
done

# Phase 4: P2跨疾病
info "═══════════════════════════════════════════════════════════════"
info "  Phase 4/5: P2跨疾病验证"
info "═══════════════════════════════════════════════════════════════"
for gse in GSE138425 GSE162429 GSE141512 GSE213677; do
    AVAIL=$(df -m "$PROJECT" | awk 'NR==2{print $4}')
    [ "$AVAIL" -lt 2048 ] && { warn "磁盘不足"; break; }
    if download_series_matrix "$gse"; then SUCCESS=$((SUCCESS+1)); else FAILED=$((FAILED+1)); fi
    download_soft "$gse"
    sleep 2
done

# Phase 5: P2机制
info "═══════════════════════════════════════════════════════════════"
info "  Phase 5/5: P2机制深化"
info "═══════════════════════════════════════════════════════════════"
for gse in GSE262714 GSE40066 GSE242046 GSE246410 GSE196192 GSE198687; do
    AVAIL=$(df -m "$PROJECT" | awk 'NR==2{print $4}')
    [ "$AVAIL" -lt 1536 ] && { warn "磁盘不足，停止"; break; }
    if download_series_matrix "$gse"; then SUCCESS=$((SUCCESS+1)); else FAILED=$((FAILED+1)); fi
    download_soft "$gse"
    sleep 2
done

# 完成
END=$(date +%s)
DUR=$((END - START))
H=$((DUR/3600)); M=$(((DUR%3600)/60))
echo ""
echo -e "${G}╔══════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${G}║     ✅ GEO批量下载完成! 成功:$SUCCESS | 失败:$FAILED                         ║${NC}"
echo -e "${G}║     耗时: ${H}小时${M}分钟                                                      ║${NC}"
echo -e "${G}╚══════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "📊 磁盘使用:"; df -h "$PROJECT"
echo ""
echo "📁 已下载数据:"; du -sh "$GEO"/* 2>/dev/null | sort -rh | head -30
echo ""
echo "⚠️ 需手动下载:"
echo "  GTEx v8: https://gtexportal.org → Heart_LV.v8.allpairs.txt.gz → $MR/"
echo "  FinnGen: https://storage.finngen.fi → finngen_R10_I9_HEARTFAIL.gz → $MR/"
