#!/bin/bash
SMR="~/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/smr_software/smr"
OUT="~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables"
DATA="~/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data"

echo "=== SMR eQTLGen → HERMES ==="
$SMR --bfile $DATA/eQTLGen_ref/1kg_eur_chr19 --gwas-summary $DATA/HERMES_HF.ma --eqtl-summary $DATA/eQTLGen_NDUFB7.ma --out $OUT/SMR_eQTLGen_HF --smr-multi --heidi --heidi-mtd 0 --thread-num 4 2>/dev/null || echo "eQTLGen SMR失败"

echo "=== SMR GTEx → HERMES ==="
$SMR --bfile $DATA/GTEx_ref/gtex_chr19 --gwas-summary $DATA/HERMES_HF.ma --eqtl-summary $DATA/GTEx_Heart_LV_NDUFB7.ma --out $OUT/SMR_GTEx_HF --smr-multi --heidi --heidi-mtd 0 --thread-num 4 2>/dev/null || echo "GTEx SMR失败"

echo "=== SMR完成 ==="
