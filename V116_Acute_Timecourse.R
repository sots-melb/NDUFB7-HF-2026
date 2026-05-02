#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mclust)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V116_Acute")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V116: 急性损伤时间动态双峰验证")
message("========================================")

test_bimodal <- function(expr_vec, label) {
  expr_vec <- as.numeric(expr_vec); expr_vec <- expr_vec[!is.na(expr_vec)]
  n <- length(expr_vec)
  if (n < 20) return(data.frame(label=label, n=n, g=NA, peak1=NA, peak2=NA, gap=NA, is_bimodal=NA, zero_pct=NA, mid_pct=NA, all_or_none=NA, stringsAsFactors=FALSE))
  mod <- tryCatch(densityMclust(expr_vec, G=1:3, verbose=FALSE), error=function(e) NULL)
  if (is.null(mod)) return(data.frame(label=label, n=n, g=NA, peak1=NA, peak2=NA, gap=NA, is_bimodal=FALSE, zero_pct=mean(expr_vec==0)*100, mid_pct=NA, all_or_none=NA, stringsAsFactors=FALSE))
  best_g <- mod$G; peaks <- if(best_g>=2) sort(mod$parameters$mean) else NA
  gap <- if(best_g>=2 && length(peaks)>=2) abs(peaks[2]-peaks[1]) else 0
  is_bimodal <- best_g==2 && gap > 0.5*sd(expr_vec)
  zero_pct <- mean(expr_vec==0)*100; nonzero <- expr_vec[expr_vec>0]
  high_thresh <- ifelse(length(nonzero)>0, quantile(nonzero,0.9), 0)
  mid_pct <- mean(expr_vec>0 & expr_vec<high_thresh)*100
  all_or_none <- zero_pct + mean(expr_vec>=high_thresh)*100
  data.frame(label=label, n=n, g=best_g, peak1=ifelse(length(peaks)>=1,peaks[1],NA), peak2=ifelse(length(peaks)>=2,peaks[2],NA), gap=gap, is_bimodal=is_bimodal, zero_pct=zero_pct, mid_pct=mid_pct, all_or_none=all_or_none, stringsAsFactors=FALSE)
}

res_list <- list()

# --- A. GSE157282 DOX toxicity (check Downloads) ---
# RAW.tar需解压，此处仅检查并提示
f157282 <- "~/Downloads/GSE157282_RAW.tar"
if (file.exists(f157282)) {
  message("[INFO] GSE157282 RAW.tar found (", round(file.size(f157282)/1e6,1), "MB). Need manual unpack for full analysis.")
  # 尝试寻找已解压的series matrix
  sm <- list.files("~/Downloads", pattern="GSE157282.*series_matrix", full.names=TRUE)[1]
  if (!is.na(sm)) {
    message("[LOAD] GSE157282 series matrix for phenotype...")
    # 可解析分组，但表达矩阵需RAW解压
  }
} else { message("[SKIP] GSE157282 RAW.tar not found") }

# --- B. GSE104150 time series ---
f104150 <- "~/Downloads/GSE104150_RAW.tar"
if (file.exists(f104150)) {
  message("[INFO] GSE104150 RAW.tar found. Checking for processed data...")
}
# 尝试读取project中可能存在的处理文件
f104_proc <- list.files("01_data", pattern="GSE104150", recursive=TRUE, full.names=TRUE)[1]
if (!is.na(f104_proc) && file.exists(f104_proc)) {
  message("[LOAD] GSE104150 processed: ", basename(f104_proc))
  # 根据实际格式调整
} else { message("[SKIP] GSE104150 processed data not found") }

# --- C. GSE214611 STEMI Visium spatial ---
# 使用project中的Visium h5ad（如果存在）或spatial CSV
f214611 <- list.files("01_data/02_spatial", pattern="GSE214611|STEMI", full.names=TRUE, recursive=TRUE)[1]
if (is.na(f214611)) f214611 <- list.files("03_results", pattern="GSE214611", full.names=TRUE, recursive=TRUE)[1]
if (!is.na(f214611) && file.exists(f214611)) {
  message("[LOAD] GSE214611 spatial: ", basename(f214611))
  # 如果是csv
  if (grepl("\\.csv$", f214611)) {
    dt <- fread(f214611)
    if ("NDUFB7" %in% colnames(dt)) {
      expr <- dt$NDUFB7
      res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE214611_STEMI_Visium")
    }
  }
} else { message("[SKIP] GSE214611 spatial not found") }

# --- D. GSE315590 mouse TAC (if in project) ---
f315590 <- list.files("01_data", pattern="GSE315590", recursive=TRUE, full.names=TRUE)[1]
if (!is.na(f315590) && file.exists(f315590)) {
  message("[LOAD] GSE315590 mouse TAC")
  obj <- readRDS(f315590)
  if (inherits(obj, "Seurat")) {
    expr <- as.numeric(FetchData(obj, vars="Ndufb7")$Ndufb7)
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE315590_Mouse_TAC")
  }
} else { message("[SKIP] GSE315590 not found") }

# --- 保存 ---
if (length(res_list)>0) {
  df_res <- do.call(rbind, res_list)
  write.csv(df_res, file.path(outdir,"V116_acute_bimodal.csv"), row.names=FALSE)
  message("\n=== V116 急性损伤汇总 ===")
  print(df_res[,c("label","n","g","is_bimodal","zero_pct","all_or_none")])
} else {
  message("\n[WARN] V116: 无可用急性损伤数据集完成分析")
  message("[ACTION] 请手动解压 GSE157282/GSE104150 RAW.tar 后重新执行")
}

message("[DONE] V116: ", outdir)
