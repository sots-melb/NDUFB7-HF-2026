#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mclust)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V117_Independent")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V117: 独立队列 + 胚胎发育双峰验证")
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

# --- A. GSE168742 independent human sc ---
# 优先使用project中已处理数据
f168 <- list.files("01_data", pattern="GSE168742.*\\.rds$|GSE168742.*\\.Robj$|168742.*\\.rds$", recursive=TRUE, full.names=TRUE)[1]
if (is.na(f168)) f168 <- list.files("03_results", pattern="GSE168742", recursive=TRUE, full.names=TRUE)[1]
if (!is.na(f168) && file.exists(f168)) {
  message("[LOAD] GSE168742: ", basename(f168))
  obj <- readRDS(f168)
  if (inherits(obj, "Seurat")) {
    expr <- as.numeric(FetchData(obj, vars="NDUFB7")$NDUFB7)
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE168742_Human_sc")
    message("[PASS] GSE168742: n=", length(expr))
  } else { message("[WARN] GSE168742 object not Seurat") }
} else { message("[SKIP] GSE168742 processed file not found") }

# --- B. GSE106118 embryo development trajectory ---
f106118 <- "~/Downloads/GSE106118_UMI_count_merge.txt.gz"
if (file.exists(f106118)) {
  message("[LOAD] GSE106118 embryo UMI...")
  dt <- fread(f106118, header=TRUE, nrows=50000)
  gene_col <- colnames(dt)[1]
  ndufb7_idx <- which(toupper(dt[[gene_col]])=="NDUFB7")
  if (length(ndufb7_idx)>0) {
    expr <- as.numeric(dt[ndufb7_idx, -1, with=FALSE])
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE106118_Embryo")
    message("[PASS] GSE106118 embryo: n=", length(expr))
  } else { message("[WARN] NDUFB7 not found in GSE106118") }
} else { message("[SKIP] GSE106118 not found") }

# --- C. GSE116250 RPKM (Bulk, if available in project) ---
f116250 <- list.files("01_data", pattern="GSE116250|116250", recursive=TRUE, full.names=TRUE)[1]
if (is.na(f116250)) f116250 <- list.files("03_results", pattern="GSE116250", recursive=TRUE, full.names=TRUE)[1]
if (!is.na(f116250) && file.exists(f116250)) {
  message("[LOAD] GSE116250: ", basename(f116250))
  # 根据实际格式处理（此处为框架）
  message("[INFO] GSE116250 loaded, format-specific extraction needed")
} else { message("[SKIP] GSE116250 not found") }

# --- D. GSE141910 (DCM bulk, if in project) ---
f141910 <- list.files("01_data", pattern="GSE141910", recursive=TRUE, full.names=TRUE)[1]
if (!is.na(f141910) && file.exists(f141910)) {
  message("[LOAD] GSE141910: ", basename(f141910))
} else { message("[SKIP] GSE141910 not found") }

# --- 保存 ---
if (length(res_list)>0) {
  df_res <- do.call(rbind, res_list)
  write.csv(df_res, file.path(outdir,"V117_independent_bimodal.csv"), row.names=FALSE)
  message("\n=== V117 独立/胚胎汇总 ===")
  print(df_res[,c("label","n","g","is_bimodal","zero_pct","mid_pct","all_or_none")])
  
  # 胚胎vs成体对比图
  p <- ggplot(df_res, aes(x=label, y=zero_pct, fill=all_or_none>70)) +
    geom_bar(stat="identity") + coord_flip() +
    scale_fill_manual(values=c("TRUE"="#FDE725","FALSE"="#31688E")) +
    labs(title="Zero%: Embryo vs Adult", subtitle="Yellow = All-or-none >70%") + theme_minimal()
  ggsave(file.path(outdir,"V117_zero_pct_comparison.png"), p, width=7, height=4, dpi=300)
} else { message("[WARN] V117: 无可用数据完成分析") }

message("[DONE] V117: ", outdir)
