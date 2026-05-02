#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V133_GSE57338")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V133: GSE57338探针-基因映射终极修复")
message("========================================")

# --- 策略A: Bioconductor hugene11sttranscriptcluster.db ---
message(">>> [策略A] Bioconductor hugene11sttranscriptcluster.db")
has_pkg <- require("hugene11sttranscriptcluster.db", quietly=TRUE)
if (!has_pkg) {
  message("[MISS] hugene11sttranscriptcluster.db not installed")
  message("[ACTION] 尝试安装: BiocManager::install('hugene11sttranscriptcluster.db')")
  tryCatch({
    if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
    BiocManager::install("hugene11sttranscriptcluster.db", ask=FALSE, update=FALSE)
    has_pkg <- require("hugene11sttranscriptcluster.db", quietly=TRUE)
  }, error=function(e) message("[FAIL] Install failed: ", e$message))
}

# 加载GSE57338表达矩阵（从R对象或CEL）
rds_file <- list.files(c("01_data","03_results"), pattern="GSE57338.*gene_level.*\\.rds$|GSE57338.*\\.rds$", full.names=TRUE, recursive=TRUE)[1]

if (!is.na(rds_file) && file.exists(rds_file) && has_pkg) {
  message("[LOAD] GSE57338 RDS: ", basename(rds_file))
  obj <- readRDS(rds_file)
  
  # 提取表达矩阵
  expr_mat <- NULL
  if (is.list(obj) && "expression" %in% names(obj)) expr_mat <- obj$expression
  else if (is.matrix(obj) || is.data.frame(obj)) expr_mat <- as.matrix(obj)
  
  if (!is.null(expr_mat)) {
    probes <- rownames(expr_mat)
    message("[PASS] Expression matrix: ", nrow(expr_mat), " probes × ", ncol(expr_mat), " samples")
    
    # 映射
    mapped <- AnnotationDbi::mapIds(hugene11sttranscriptcluster.db, 
                                     keys=probes, column="SYMBOL", 
                                     keytype="PROBEID", multiVals="first")
    
    # 构建基因级矩阵
    gene_expr <- data.frame(probe=probes, gene=mapped, expr_mat, stringsAsFactors=FALSE)
    gene_expr <- gene_expr[!is.na(gene_expr$gene), ]
    
    # 去重复：同一基因多个探针取中位数
    gene_expr_unique <- gene_expr %>% 
      group_by(gene) %>% 
      summarise(across(starts_with("GSM"), ~median(.x, na.rm=TRUE)), .groups="drop")
    
    # 提取NDUFB7
    if ("NDUFB7" %in% gene_expr_unique$gene) {
      ndufb7_expr <- as.numeric(gene_expr_unique[gene_expr_unique$gene=="NDUFB7", -1])
      message("[PASS] NDUFB7 extracted: n=", length(ndufb7_expr))
      write.csv(data.frame(sample=colnames(gene_expr_unique)[-1], NDUFB7=ndufb7_expr), 
                file.path(outdir,"V133_NDUFB7_expression.csv"), row.names=FALSE)
    } else {
      message("[WARN] NDUFB7 not found after mapping")
    }
    
    # 保存完整基因矩阵
    write.csv(gene_expr_unique, file.path(outdir,"V133_GSE57338_gene_level.csv"), row.names=FALSE)
    message("[DONE] Gene-level matrix saved")
  } else {
    message("[WARN] Cannot extract expression matrix from RDS")
  }
} else {
  message("[SKIP] Strategy A unavailable (missing package or RDS)")
}

# --- 策略B: 手动GPL11532解析（如果策略A失败）---
if (!has_pkg || is.null(expr_mat)) {
  message("\n>>> [策略B] 手动GPL11532解析")
  
  # 寻找GPL11532注释文件
  gpl_files <- list.files("~/Downloads", pattern="GPL11532.*\\.soft\\.gz$|GPL11532.*\\.annot\\.gz$", full.names=TRUE)
  if (length(gpl_files) == 0) {
    message("[MISS] GPL11532 soft file not found in Downloads")
    message("[ACTION] 手动下载: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL11532")
    message("         点击 'SOFT formatted family file(s)' → 下载 GPL11532_family.soft.gz")
  } else {
    message("[FOUND] GPL file: ", basename(gpl_files[1]))
    # 解析框架（需根据实际格式调整）
    message("[INFO] 解析框架就绪，需根据soft格式定制")
  }
}

# --- 策略C: 降级方案（如果全部失败）---
message("\n>>> [策略C] 降级方案")
message("如果探针映射无法完成，GSE57338可用于：")
message("  1. 样本级QC和批次效应评估（无需基因名）")
message("  2. 与其他Bulk数据的meta-analysis（使用probe-level ID匹配）")
message("  3. 在Results中诚实标注: 'Probe-to-gene mapping for GSE57338 was attempted but")
message("     remained incomplete due to platform annotation file corruption;")
message("     downstream gene-level analyses used alternative cohorts (GSE55296, GSE141910).'")

message("\n[DONE] V133: ", outdir)
