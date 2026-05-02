#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mclust)
  library(moments)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V115_Cross_Model")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V115: 跨病因/跨物种双峰验证")
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

# --- A. GSE55296 Bulk ICM/DCM/Control ---
f55296 <- "~/Downloads/GSE55296_count_data.txt.gz"
if (file.exists(f55296)) {
  message("[LOAD] GSE55296 count data...")
  dt <- fread(f55296, header=TRUE)
  # 假设第一列是基因名
  gene_col <- colnames(dt)[1]
  ndufb7_row <- dt[[gene_col]][toupper(dt[[gene_col]])=="NDUFB7"]
  if (length(ndufb7_row)>0) {
    expr <- as.numeric(dt[toupper(dt[[gene_col]])=="NDUFB7", -1, with=FALSE])
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE55296_Bulk_All")
    message("[PASS] GSE55296 NDUFB7: n=", length(expr), " mean=", round(mean(expr),2))
  } else { message("[WARN] NDUFB7 not found in GSE55296 rownames") }
} else { message("[SKIP] GSE55296 not in Downloads") }

# --- B. GSE121893 human heart sc ---
f121893 <- "~/Downloads/GSE121893_human_heart_sc_umi.csv.gz"
if (file.exists(f121893)) {
  message("[LOAD] GSE121893 sc UMI (sampling if large)...")
  dt <- fread(f121893, header=TRUE, nrows=50000)  # 先读前50k行探测结构
  # 假设行为基因，列为细胞；或相反。取第一列作为基因ID
  gene_col <- colnames(dt)[1]
  ndufb7_idx <- which(toupper(dt[[gene_col]])=="NDUFB7")
  if (length(ndufb7_idx)>0) {
    expr <- as.numeric(dt[ndufb7_idx, -1, with=FALSE])
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE121893_sc_All")
    message("[PASS] GSE121893 NDUFB7: n=", length(expr))
  } else { message("[WARN] NDUFB7 not found in GSE121893") }
} else { message("[SKIP] GSE121893 not found") }

# --- C. GSE275031 mouse HFpEF ---
f275031 <- "~/Downloads/GSE275031_integrated_seurat_obj.rds.gz"
if (file.exists(f275031)) {
  message("[LOAD] GSE275031 mouse HFpEF Seurat...")
  srt <- readRDS(f275031)
  # 小鼠基因名: Ndufb7
  if ("Ndufb7" %in% rownames(srt)) {
    expr <- as.numeric(FetchData(srt, vars="Ndufb7")$Ndufb7)
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE275031_Mouse_HFpEF")
    message("[PASS] Mouse Ndufb7: n=", length(expr))
  } else if ("NDUFB7" %in% rownames(srt)) {
    expr <- as.numeric(FetchData(srt, vars="NDUFB7")$NDUFB7)
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE275031_Mouse_HFpEF_humanID")
    message("[PASS] Mouse NDUFB7 (human ID): n=", length(expr))
  } else { message("[WARN] Ndufb7/NDUFB7 not in GSE275031") }
} else { message("[SKIP] GSE275031 not found") }

# --- D. GSE168742 human sc (if processed in project) ---
f168742 <- list.files("01_data", pattern="GSE168742.*\\.rds$|GSE168742.*\\.Robj$", recursive=TRUE, full.names=TRUE)[1]
if (is.na(f168742)) f168742 <- list.files("03_results", pattern="GSE168742", recursive=TRUE, full.names=TRUE)[1]
if (!is.na(f168742) && file.exists(f168742)) {
  message("[LOAD] GSE168742 from project...")
  obj <- readRDS(f168742)
  # 尝试多种对象类型
  if (inherits(obj, "Seurat")) {
    expr <- as.numeric(FetchData(obj, vars="NDUFB7")$NDUFB7)
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE168742_sc")
  } else { message("[WARN] GSE168742 object type unknown") }
} else { message("[SKIP] GSE168742 processed file not found") }

# --- 保存 ---
if (length(res_list)>0) {
  df_res <- do.call(rbind, res_list)
  write.csv(df_res, file.path(outdir,"V115_cross_model_bimodal.csv"), row.names=FALSE)
  message("\n=== V115 跨模型双峰汇总 ===")
  print(df_res[,c("label","n","g","is_bimodal","zero_pct","mid_pct","all_or_none")])
  
  # 跨模型对比图
  p <- ggplot(df_res, aes(x=label, y=all_or_none, fill=is_bimodal)) + 
    geom_bar(stat="identity") + coord_flip() + 
    scale_fill_manual(values=c("TRUE"="#FDE725","FALSE"="#440154","NA"="grey50")) +
    labs(title="All-or-None Index Across Models", x="", y="Zero% + Top10%") + theme_minimal()
  ggsave(file.path(outdir,"V115_cross_model_bar.png"), p, width=8, height=5, dpi=300)
} else { message("[WARN] No successful bimodal tests") }

message("[DONE] V115: ", outdir)
