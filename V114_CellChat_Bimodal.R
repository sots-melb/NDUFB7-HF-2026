#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(mclust)
  library(moments)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V114_Bimodal")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V114: CellChat提取 + NDUFB7双峰检验")
message("========================================")

test_bimodal <- function(expr_vec, label) {
  expr_vec <- as.numeric(expr_vec); expr_vec <- expr_vec[!is.na(expr_vec)]
  n <- length(expr_vec)
  if (n < 30) return(data.frame(label=label, n=n, g=NA, peak1=NA, peak2=NA, gap=NA, sd=NA, is_bimodal=NA, zero_pct=NA, mid_pct=NA, all_or_none=NA, stringsAsFactors=FALSE))
  mod <- tryCatch(densityMclust(expr_vec, G=1:3, verbose=FALSE), error=function(e) NULL)
  if (is.null(mod)) return(data.frame(label=label, n=n, g=NA, peak1=NA, peak2=NA, gap=NA, sd=sd(expr_vec), is_bimodal=FALSE, zero_pct=mean(expr_vec==0)*100, mid_pct=NA, all_or_none=NA, stringsAsFactors=FALSE))
  best_g <- mod$G; peaks <- if(best_g>=2) sort(mod$parameters$mean) else NA
  gap <- if(best_g>=2 && length(peaks)>=2) abs(peaks[2]-peaks[1]) else 0
  sd_expr <- sd(expr_vec); is_bimodal <- best_g==2 && gap > 0.5*sd_expr
  zero_pct <- mean(expr_vec==0)*100; nonzero <- expr_vec[expr_vec>0]
  high_thresh <- ifelse(length(nonzero)>0, quantile(nonzero,0.9), 0)
  mid_pct <- mean(expr_vec>0 & expr_vec<high_thresh)*100
  all_or_none <- zero_pct + mean(expr_vec>=high_thresh)*100
  data.frame(label=label, n=n, g=best_g, peak1=ifelse(length(peaks)>=1,peaks[1],NA), peak2=ifelse(length(peaks)>=2,peaks[2],NA), gap=gap, sd=sd_expr, is_bimodal=is_bimodal, zero_pct=zero_pct, mid_pct=mid_pct, all_or_none=all_or_none, stringsAsFactors=FALSE)
}

res_list <- list()
# --- GSE183852 CM ---
cm_file <- list.files("01_data", pattern="CM.*annotated.*\\.rds$|Pure_CM.*\\.RDS$|34_srt.*\\.rds$", recursive=TRUE, full.names=TRUE)[1]
if (!is.na(cm_file) && file.exists(cm_file)) {
  srt <- readRDS(cm_file); message("[PASS] CM: ", ncol(srt), " cells from ", basename(cm_file))
  ndufb7 <- as.numeric(FetchData(srt, vars="NDUFB7")$NDUFB7)
  res_list[[1]] <- test_bimodal(ndufb7, "GSE183852_CM_All")
  if ("condition" %in% colnames(srt@meta.data)) {
    for (cond in unique(srt$condition)) {
      cells <- colnames(srt)[srt$condition==cond]
      if (length(cells)>30) res_list[[length(res_list)+1]] <- test_bimodal(ndufb7[cells], paste0("GSE183852_CM_",cond))
    }
  }
  if ("NDUFB7_level" %in% colnames(srt@meta.data)) {
    for (lvl in unique(srt$NDUFB7_level)) {
      cells <- colnames(srt)[srt$NDUFB7_level==lvl]
      if (length(cells)>10) res_list[[length(res_list)+1]] <- test_bimodal(ndufb7[cells], paste0("GSE183852_",lvl))
    }
  }
  p <- ggplot(data.frame(NDUFB7=ndufb7), aes(x=NDUFB7)) + geom_density(fill="#31688E", alpha=0.3) + geom_rug(alpha=0.1) + labs(title="NDUFB7 Distribution in CM", subtitle=paste0("All-or-none: ", round(res_list[[1]]$all_or_none,1), "%")) + theme_minimal()
  ggsave(file.path(outdir,"V114_GSE183852_density.png"), p, width=6, height=4, dpi=300)
} else { message("[FAIL] No CM file found") }

# --- CellChat ---
chat_dir <- "03_results/V113A_CellChat"
if (dir.exists(chat_dir)) {
  f <- list.files(chat_dir, pattern="all_pairs\\.csv$", full.names=TRUE)[1]
  if (!is.na(f)) {
    df <- read.csv(f)
    message("\n=== CellChat Top 5 ==="); print(head(df[,c("source","target","ligand","receptor","prob")], 5))
    write.csv(df, file.path(outdir,"V114_cellchat_matrix.csv"), row.names=FALSE)
  }
} else { message("[SKIP] V113A not found") }

if (length(res_list)>0) {
  df_res <- do.call(rbind, res_list)
  write.csv(df_res, file.path(outdir,"V114_bimodal_summary.csv"), row.names=FALSE)
  message("\n=== V114 双峰汇总 ==="); print(df_res[,c("label","n","g","is_bimodal","zero_pct","mid_pct","all_or_none")])
}
message("[DONE] V114: ", outdir)
