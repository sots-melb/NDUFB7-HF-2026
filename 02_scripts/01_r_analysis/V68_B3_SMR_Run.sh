#!/bin/bash
SMR="$HOME/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/smr_software/smr"
OUT="$HOME/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables"
DATA="$HOME/Projects/NDUFB7_HF_2026_04_20/01_data/03_mr_data/smr_input"

echo "=== SMR eQTLGen → HERMES ==="
$SMR --bfile $DATA/1kg_eur_chr19 --gwas-summary $DATA/HERMES_HF_SMR.ma --eqtl-summary $DATA/eQTLGen_NDUFB7.ma --out $OUT/SMR_eQTLGen_HF --smr-multi --heidi --heidi-mtd 0 --thread-num 4 2>&1 | tail -20

echo "=== SMR GTEx → HERMES ==="
$SMR --bfile $DATA/1kg_eur_chr19 --gwas-summary $DATA/HERMES_HF_SMR.ma --eqtl-summary $DATA/GTEx_HeartLV_NDUFB7.ma --out $OUT/SMR_GTEx_HF --smr-multi --heidi --heidi-mtd 0 --thread-num 4 2>&1 | tail -20

echo "=== SMR完成 ==="
ls -lh $OUT/SMR_*
