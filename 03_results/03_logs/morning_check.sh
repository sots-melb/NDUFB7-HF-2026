#!/bin/bash
PROJECT="$HOME/Projects/NDUFB7_HF_2026_04_20"
NIGHT_DIR=$(ls -td "$PROJECT/03_results/03_logs/night_batch_"* 2>/dev/null | head -1)
LOG="$NIGHT_DIR/MASTER.log"

echo "========================================"
echo "V62 明早检查 | $(date)"
echo "========================================"

echo ""
echo "[1] 夜间任务总状态"
grep -E "SUCCESS|FAILED|GAVE UP" "$LOG" 2>/dev/null | tail -20

echo ""
echo "[2] 下载文件审计"
for d in GSE315590 GSE59867; do
    dir="$PROJECT/01_data/01_raw_geo/$d"
    echo "  $d:"
    ls -lh "$dir" 2>/dev/null | awk '{print "    " $5, $9}' || echo "    无文件"
done

echo "  Visium:"
ls -lh "$PROJECT/01_data/02_spatial/"*.h5ad 2>/dev/null | awk '{print "    " $5, $9}' || echo "    无h5ad"

echo "  GPL11532:"
ls -lh "$PROJECT/01_data/01_raw_geo/GSE57338/GPL11532_family.soft.gz" 2>/dev/null | awk '{print "    " $5, $9}' || echo "    无文件"

echo ""
echo "[3] R包安装状态"
Rscript -e 'pkgs <- c("spdep","Seurat","Scissor","spacexr","scTenifoldNet"); for(p in pkgs) cat(p, ":", require(p, character.only=TRUE, quietly=TRUE), "\n")'

echo ""
echo "[4] Python包安装状态"
python3 -c "
import importlib
for pkg in ['scFEA', 'cell2location', 'scanpy']:
    try:
        m = importlib.import_module(pkg)
        print(f'  {pkg}: OK')
    except:
        print(f'  {pkg}: FAIL')
"

echo ""
echo "[5] SMR二进制"
ls -lh "$PROJECT/01_data/03_mr_data/smr_software/smr" 2>/dev/null || echo "  SMR: 未找到"

echo ""
echo "[6] 活跃进程"
ps aux | grep -E "wget|Rscript|nohup" | grep -v grep | wc -l | xargs -I{} echo "  后台进程数: {}"

echo ""
echo "========================================"
echo "检查完成: $(date)"
echo "========================================"
