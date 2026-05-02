#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(mclust)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V130_GSE168742")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V130: GSE168742独立队列验证")
message("========================================")

# 加载GSE168742（之前V117已确认存在）
f168 <- list.files("01_data", pattern="GSE168742.*\\.rds$|168742.*\\.rds$", recursive=TRUE, full.names=TRUE)[1]
if (is.na(f168)) f168 <- list.files("03_results", pattern="GSE168742", recursive=TRUE, full.names=TRUE)[1]
if (is.na(f168) || !file.exists(f168)) {
  message("[FAIL] GSE168742 not found"); quit(status=1)
}

message("[LOAD] ", basename(f168))
obj <- readRDS(f168)

# 暴力提取NDUFB7
expr_168 <- NULL
try_list <- list(
  tryCatch(as.numeric(obj@assays$RNA@counts["NDUFB7", ]), error=function(e) NULL),
  tryCatch(as.numeric(obj@assays$RNA$counts["NDUFB7", ]), error=function(e) NULL),
  tryCatch(as.numeric(obj@assays$RNA@data["NDUFB7", ]), error=function(e) NULL),
  tryCatch(as.numeric(obj[["RNA"]]$counts["NDUFB7", ]), error=function(e) NULL)
)
for (i in seq_along(try_list)) {
  if (!is.null(try_list[[i]]) && length(try_list[[i]]) > 10) {
    expr_168 <- try_list[[i]]
    message("[FOUND] Method ", i, ": n=", length(expr_168))
    break
  }
}
# 如果都不行，尝试矩阵行列探测
if (is.null(expr_168) && (is.matrix(obj) || is.data.frame(obj))) {
  rnames <- rownames(obj)
  if (!is.null(rnames) && "NDUFB7" %in% rnames) {
    expr_168 <- as.numeric(obj["NDUFB7", ])
    message("[FOUND] Matrix row: n=", length(expr_168))
  } else {
    cnames <- colnames(obj)
    if (!is.null(cnames) && "NDUFB7" %in% cnames) {
      expr_168 <- as.numeric(obj[, "NDUFB7"])
      message("[FOUND] Matrix col: n=", length(expr_168))
    }
  }
}
if (is.null(expr_168)) {
  message("[FAIL] Cannot extract NDUFB7. Object class: ", paste(class(obj), collapse="/"))
  quit(status=1)
}

expr_168 <- expr_168[!is.na(expr_168)]
message("[PASS] GSE168742: n=", length(expr_168), " zero_pct=", round(mean(expr_168==0)*100,1), "%")

# 双峰检验
mod <- tryCatch(densityMclust(expr_168, G=1:3, verbose=FALSE), error=function(e) NULL)
if (!is.null(mod)) {
  best_g <- mod$G
  peaks <- if(best_g>=2) sort(mod$parameters$mean) else NA
  zero_pct <- mean(expr_168==0)*100
  nonzero <- expr_168[expr_168>0]
  high_thresh <- ifelse(length(nonzero)>0, quantile(nonzero,0.9), 0)
  all_or_none <- zero_pct + mean(expr_168>=high_thresh)*100
  
  res <- data.frame(
    Dataset="GSE168742", n=length(expr_168), G=best_g, 
    zero_pct=zero_pct, all_or_none=all_or_none,
    peak1=ifelse(length(peaks)>=1,peaks[1],NA),
    peak2=ifelse(length(peaks)>=2,peaks[2],NA)
  )
  write.csv(res, file.path(outdir, "V168742_bimodal_stats.csv"), row.names=FALSE)
  message("\n=== GSE168742 分布 ===")
  print(res)
  
  # 密度图
  p <- ggplot(data.frame(NDUFB7=expr_168), aes(x=NDUFB7)) +
    geom_density(fill="#31688E", alpha=0.3, color="#440154", linewidth=1) +
    labs(title="GSE168742 NDUFB7 Distribution", subtitle=paste0("G=",best_g," | All-or-none=",round(all_or_none,1),"%")) +
    theme_minimal()
  ggsave(file.path(outdir,"V130_GSE168742_density.png"), p, width=6, height=4, dpi=300)
  
  # 判定
  if (best_g == 3) {
    message("[PASS] GSE168742 also shows tri-modal distribution — CONSISTENT with GSE183852")
    message("[NARRATIVE] 'Tri-modality is conserved across independent human heart failure cohorts'")
  } else if (best_g == 2) {
    message("[PARTIAL] GSE168742 is bi-modal — supports 'platform/etiology-dependent modality'")
  } else {
    message("[INFO] GSE168742 is unimodal — may reflect different disease stage or sampling")
  }
} else {
  message("[WARN] densityMclust failed on GSE168742")
}

message("\n[DONE] V130: ", outdir)
