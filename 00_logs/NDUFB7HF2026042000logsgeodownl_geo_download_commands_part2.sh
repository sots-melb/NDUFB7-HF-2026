#!/bin/bash
set -e
PROJECT='/home/y411869/Projects/NDUFB7_HF_2026_04_20'
GEO='/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo'
LOG='/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs'
mkdir -p '$GEO' '$LOG'

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE277721'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE277721/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE277721/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar' ]; then
  echo '[SKIP] GSE277721 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar'
  echo 'SKIP GSE277721 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE277721 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE277721/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE277721.log' 2>&1; then
    echo '[OK] GSE277721 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE277721/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE277721 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE277721 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar'
    echo 'FAILED GSE277721 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE277721/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277721/suppl/GSE277721_RAW.tar'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE190132'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE190132/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE190132/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar' ]; then
  echo '[SKIP] GSE190132 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar'
  echo 'SKIP GSE190132 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE190132 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE190132/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE190132.log' 2>&1; then
    echo '[OK] GSE190132 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE190132/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE190132 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE190132 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar'
    echo 'FAILED GSE190132 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE190132/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190132/suppl/GSE190132_RAW.tar'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198699'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198699/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198699/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz' ]; then
  echo '[SKIP] GSE198699 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz'
  echo 'SKIP GSE198699 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE198699 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198699/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE198699.log' 2>&1; then
    echo '[OK] GSE198699 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198699/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz | cut -f1)'
    echo 'SUCCESS GSE198699 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE198699 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz'
    echo 'FAILED GSE198699 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198699/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198699/suppl/GSE198699_gene_count.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz' ]; then
  echo '[SKIP] GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz'
  echo 'SKIP GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE184494.log' 2>&1; then
    echo '[OK] GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz | cut -f1)'
    echo 'SUCCESS GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz'
    echo 'FAILED GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_Non-normalized_data.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar' ]; then
  echo '[SKIP] GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar'
  echo 'SKIP GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE184494.log' 2>&1; then
    echo '[OK] GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar'
    echo 'FAILED GSE184494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE184494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184494/suppl/GSE184494_RAW.tar'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE196192'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE196192/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE196192/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz' ]; then
  echo '[SKIP] GSE196192 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz'
  echo 'SKIP GSE196192 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE196192 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE196192/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE196192.log' 2>&1; then
    echo '[OK] GSE196192 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE196192/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz | cut -f1)'
    echo 'SUCCESS GSE196192 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE196192 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz'
    echo 'FAILED GSE196192 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE196192/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196192/suppl/GSE196192_Lyme-PBMC-transcriptome-set2.counts.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198687'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198687/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198687/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar' ]; then
  echo '[SKIP] GSE198687 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar'
  echo 'SKIP GSE198687 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE198687 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198687/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE198687.log' 2>&1; then
    echo '[OK] GSE198687 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198687/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE198687 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE198687 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar'
    echo 'FAILED GSE198687 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE198687/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198687/suppl/GSE198687_RAW.tar'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz' ]; then
  echo '[SKIP] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz'
  echo 'SKIP GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE183852.log' 2>&1; then
    echo '[OK] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz | cut -f1)'
    echo 'SUCCESS GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz'
    echo 'FAILED GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Cells.Robj.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz' ]; then
  echo '[SKIP] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz'
  echo 'SKIP GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE183852.log' 2>&1; then
    echo '[OK] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz | cut -f1)'
    echo 'SUCCESS GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz'
    echo 'FAILED GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Integrated.Robj.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz' ]; then
  echo '[SKIP] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz'
  echo 'SKIP GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE183852.log' 2>&1; then
    echo '[OK] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz | cut -f1)'
    echo 'SUCCESS GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz'
    echo 'FAILED GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_DCM_Nuclei.Robj.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz' ]; then
  echo '[SKIP] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz'
  echo 'SKIP GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE183852.log' 2>&1; then
    echo '[OK] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz | cut -f1)'
    echo 'SUCCESS GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz'
    echo 'FAILED GSE183852 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183852/suppl/GSE183852_Integrated_Counts.csv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz' ]; then
  echo '[SKIP] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz'
  echo 'SKIP GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE120064.log' 2>&1; then
    echo '[OK] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz | cut -f1)'
    echo 'SUCCESS GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz'
    echo 'FAILED GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_clean_cell_info_summary.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz' ]; then
  echo '[SKIP] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz'
  echo 'SKIP GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE120064.log' 2>&1; then
    echo '[OK] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz | cut -f1)'
    echo 'SUCCESS GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz'
    echo 'FAILED GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_TAC_raw_umi_matrix.csv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz' ]; then
  echo '[SKIP] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz'
  echo 'SKIP GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE120064.log' 2>&1; then
    echo '[OK] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz | cut -f1)'
    echo 'SUCCESS GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz'
    echo 'FAILED GSE120064 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE120064/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120064/suppl/GSE120064_files.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz' ]; then
  echo '[SKIP] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz'
  echo 'SKIP GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE165838.log' 2>&1; then
    echo '[OK] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz | cut -f1)'
    echo 'SUCCESS GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz'
    echo 'FAILED GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_barcodes.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz' ]; then
  echo '[SKIP] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz'
  echo 'SKIP GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE165838.log' 2>&1; then
    echo '[OK] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz'
    echo 'FAILED GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_counts.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz' ]; then
  echo '[SKIP] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz'
  echo 'SKIP GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE165838.log' 2>&1; then
    echo '[OK] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz'
    echo 'FAILED GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz' ]; then
  echo '[SKIP] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz'
  echo 'SKIP GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE165838.log' 2>&1; then
    echo '[OK] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz | cut -f1)'
    echo 'SUCCESS GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz'
    echo 'FAILED GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_metadata.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz' ]; then
  echo '[SKIP] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz'
  echo 'SKIP GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE165838.log' 2>&1; then
    echo '[OK] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz | cut -f1)'
    echo 'SUCCESS GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz'
    echo 'FAILED GSE165838 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE165838/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165838/suppl/GSE165838_CARE_RNA_umap.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz' ]; then
  echo '[SKIP] GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz'
  echo 'SKIP GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE46224.log' 2>&1; then
    echo '[OK] GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz | cut -f1)'
    echo 'SUCCESS GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz'
    echo 'FAILED GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_RNASeq.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz' ]; then
  echo '[SKIP] GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz'
  echo 'SKIP GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE46224.log' 2>&1; then
    echo '[OK] GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz | cut -f1)'
    echo 'SUCCESS GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz'
    echo 'FAILED GSE46224 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE46224/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46224/suppl/GSE46224_Yang_et_al_human_heart_miRNASeq.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE48166'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE48166/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE48166/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz' ]; then
  echo '[SKIP] GSE48166 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz'
  echo 'SKIP GSE48166 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE48166 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE48166/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE48166.log' 2>&1; then
    echo '[OK] GSE48166 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE48166/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz | cut -f1)'
    echo 'SUCCESS GSE48166 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE48166 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz'
    echo 'FAILED GSE48166 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE48166/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48166/suppl/GSE48166_Cufflinks_FPKM.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_cell.metadata.csv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample1.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample12.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample13.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample15.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample17.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample2.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample27.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample28.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample29.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample30.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample32.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample33.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample34.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample39.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample4.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample41.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample42.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample5.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample6.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample7.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample8.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.barcodes.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.features.tsv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz' ]; then
  echo '[SKIP] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz'
  echo 'SKIP GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE217494.log' 2>&1; then
    echo '[OK] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz | cut -f1)'
    echo 'SUCCESS GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz'
    echo 'FAILED GSE217494 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE217494/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217494/suppl/GSE217494_sample9.matrix.mtx.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar' ]; then
  echo '[SKIP] GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar'
  echo 'SKIP GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE270788.log' 2>&1; then
    echo '[OK] GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar'
    echo 'FAILED GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_RAW.tar'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz' ]; then
  echo '[SKIP] GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz'
  echo 'SKIP GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE270788.log' 2>&1; then
    echo '[OK] GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz | cut -f1)'
    echo 'SUCCESS GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz'
    echo 'FAILED GSE270788 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE270788/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270788/suppl/GSE270788_metadata.csv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE138425'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE138425/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE138425/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar' ]; then
  echo '[SKIP] GSE138425 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar'
  echo 'SKIP GSE138425 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE138425 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE138425/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE138425.log' 2>&1; then
    echo '[OK] GSE138425 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE138425/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE138425 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE138425 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar'
    echo 'FAILED GSE138425 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE138425/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138425/suppl/GSE138425_RAW.tar'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE162429'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE162429/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE162429/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar' ]; then
  echo '[SKIP] GSE162429 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar'
  echo 'SKIP GSE162429 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE162429 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE162429/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE162429.log' 2>&1; then
    echo '[OK] GSE162429 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE162429/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE162429 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE162429 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar'
    echo 'FAILED GSE162429 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE162429/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162429/suppl/GSE162429_RAW.tar'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE141512'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE141512/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE141512/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar' ]; then
  echo '[SKIP] GSE141512 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar'
  echo 'SKIP GSE141512 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE141512 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE141512/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE141512.log' 2>&1; then
    echo '[OK] GSE141512 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE141512/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE141512 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE141512 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar'
    echo 'FAILED GSE141512 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE141512/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141512/suppl/GSE141512_RAW.tar'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar' ]; then
  echo '[SKIP] GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar'
  echo 'SKIP GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE213677.log' 2>&1; then
    echo '[OK] GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar'
    echo 'FAILED GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_RAW.tar'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz' ]; then
  echo '[SKIP] GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz'
  echo 'SKIP GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE213677.log' 2>&1; then
    echo '[OK] GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz | cut -f1)'
    echo 'SUCCESS GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz'
    echo 'FAILED GSE213677 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE213677/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213677/suppl/GSE213677_gene_count_matrix.csv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE262714'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE262714/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE262714/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz' ]; then
  echo '[SKIP] GSE262714 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz'
  echo 'SKIP GSE262714 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE262714 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE262714/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE262714.log' 2>&1; then
    echo '[OK] GSE262714 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE262714/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz | cut -f1)'
    echo 'SUCCESS GSE262714 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE262714 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz'
    echo 'FAILED GSE262714 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE262714/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE262nnn/GSE262714/suppl/GSE262714_raw_counts_heart.csv.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE40066'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE40066/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE40066/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar' ]; then
  echo '[SKIP] GSE40066 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar'
  echo 'SKIP GSE40066 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE40066 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE40066/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE40066.log' 2>&1; then
    echo '[OK] GSE40066 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE40066/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE40066 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE40066 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar'
    echo 'FAILED GSE40066 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE40066/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40066/suppl/GSE40066_RAW.tar'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE242046'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE242046/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE242046/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz' ]; then
  echo '[SKIP] GSE242046 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz'
  echo 'SKIP GSE242046 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE242046 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE242046/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE242046.log' 2>&1; then
    echo '[OK] GSE242046 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE242046/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz | cut -f1)'
    echo 'SUCCESS GSE242046 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE242046 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz'
    echo 'FAILED GSE242046 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE242046/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242046/suppl/GSE242046_matrix.count.txt.gz'
  fi
fi

mkdir -p '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE246410'
if [ -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE246410/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar' ] && [ -s '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE246410/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar' ]; then
  echo '[SKIP] GSE246410 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar'
  echo 'SKIP GSE246410 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
else
  echo '[DOWN] GSE246410 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar ...'
  if wget -c -t 10 --timeout=600 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar' -O '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE246410/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/GSE246410.log' 2>&1; then
    echo '[OK] GSE246410 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar $(du -h /home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE246410/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar | cut -f1)'
    echo 'SUCCESS GSE246410 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
  else
    echo '[FAIL] GSE246410 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar'
    echo 'FAILED GSE246410 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar $(date)' >> '/home/y411869/Projects/NDUFB7_HF_2026_04_20/04_logs/master_download.log'
    rm -f '/home/y411869/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE246410/ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246410/suppl/GSE246410_RAW.tar'
  fi
fi

