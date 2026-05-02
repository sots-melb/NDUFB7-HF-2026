#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(mclust)
  library(data.table)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V132_GSE168742")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V132: GSE168742 data.table适配修复")
message("========================================")

f168 <- list.files("01_data", pattern="GSE168742.*\\.rds$|168742.*\\.rds$", recursive=TRUE, full.names=TRUE)[1]
if (is.na(f168)) f168 <- list.files("03_results", pattern="GSE168742", recursive=TRUE, full.names=TRUE)[1]
if (is.na(f168) || !file.exists(f168)) { message("[FAIL] Not found"); quit(status=1) }

message("[LOAD] ", basename(f168))
obj <- readRDS(f168)
message("[INFO] Class: ", paste(class(obj), collapse="/"))

# data.table/data.frame适配
expr_168 <- NULL
if (inherits(obj, "data.table") || is.data.frame(obj)) {
  # 检查列名是否含NDUFB7
  cnames <- colnames(obj)
  if ("NDUFB7" %in% cnames) {
    expr_168 <- as.numeric(obj[[ "NDUFB7" ]])
    message("[FOUND] Column NDUFB7: n=", length(expr_168))
  } else {
    # 检查是否有基因名在某一列
    message("[INFO] Columns: ", paste(head(cnames, 10), collapse=", "))
    # 可能第一列是基因名，其余是样本
    first_col <- as.character(obj[[1]])
    idx <- which(toupper(first_col) == "NDUFB7")[1]
    if (!is.na(idx)) {
      # 提取该行所有样本值
      expr_168 <- as.numeric(obj[idx, -1, with=FALSE])
      message("[FOUND] Row NDUFB7: n=", length(expr_168))
    }
  }
}

if (is.null(expr_168)) {
  message("[FAIL] Cannot find NDUFB7 in data.table/data.frame")
  message("[ACTION] Manual inspection needed: str(obj)")
  quit(status=1)
}

expr_168 <- expr_168[!is.na(expr_168)]
message("[PASS] Final n=", length(expr_168), " zero_pct=", round(mean(expr_168==0)*100,1), "%")

# 双峰检验
mod <- tryCatch(densityMclust(expr_168, G=1:3, verbose=FALSE), error=function(e) NULL)
if (!is.null(mod)) {
  best_g <- mod$G
  peaks <- if(best_g>=2) sort(mod$parameters$mean) else NA
  zero_pct <- mean(expr_168==0)*100
  nonzero <- expr_168[expr_168>0]
  high_thresh <- ifelse(length(nonzero)>0, quantile(nonzero,0.9), 0)
  all_or_none <- zero_pct + mean(expr_168>=high_thresh)*100
  
  res <- data.frame(Dataset="GSE168742", n=length(expr_168), G=best_g, 
                    zero_pct=zero_pct, all_or_none=all_or_none,
                    peak1=ifelse(length(peaks)>=1,peaks[1],NA),
                    peak2=ifelse(length(peaks)>=2,peaks[2],NA),
                    peak3=ifelse(length(peaks)>=3,peaks[3],NA))
  write.csv(res, file.path(outdir,"V132_GSE168742_bimodal.csv"), row.names=FALSE)
  message("\n=== GSE168742 ===")
  print(res)
  
  p <- ggplot(data.frame(NDUFB7=expr_168), aes(x=NDUFB7)) +
    geom_density(fill="#31688E", alpha=0.3, color="#440154", linewidth=1) +
    labs(title="GSE168742 NDUFB7 Distribution", subtitle=paste0("G=",best_g," | All-or-none=",round(all_or_none,1),"%")) +
    theme_minimal()
  ggsave(file.path(outdir,"V132_GSE168742_density.png"), p, width=6, height=4, dpi=300)
  
  if (best_g == 3) message("[PASS] Tri-modal — CONSISTENT with GSE183852")
  else if (best_g == 2) message("[PARTIAL] Bi-modal — platform/etiology dependent")
  else message("[INFO] Uni-modal — may reflect control-enriched sampling")
} else { message("[WARN] densityMclust failed") }

message("[DONE] V132: ", outdir)
