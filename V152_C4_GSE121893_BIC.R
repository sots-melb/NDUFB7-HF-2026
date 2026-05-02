#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(mclust) })
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V152_C4_GSE121893_BIC")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V152: C4强化 — GSE121893 BIC验证")
message("========================================")

# === 1. 搜索GSE121893 NDUFB7 ===
message("\n[1/3] 搜索GSE121893 NDUFB7...")
e121 <- NULL

# 路径1: CSV
csv_candidates <- list.files(c("03_results","01_data"), pattern = "GSE121893.*NDUFB7.*\\.csv", full.names = TRUE, recursive = TRUE)
for (f in csv_candidates) {
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (!is.null(d)) {
    nc <- names(d)[sapply(d, is.numeric)]
    if (length(nc) > 0) {
      vec <- as.numeric(unlist(d[, ..nc]))
      vec <- vec[is.finite(vec)]
      if (length(vec) > 100) { e121 <- vec; message("[PASS] CSV: ", f); break }
    }
  }
}

# 路径2: h5ad Python提取
if (is.null(e121)) {
  h5_candidates <- list.files("01_data", pattern = "GSE121893.*\\.h5ad", full.names = TRUE, recursive = TRUE)
  if (length(h5_candidates) > 0) {
    tmpcsv <- file.path(outdir, "tmp_gse121893_ndufb7.csv")
    cmd <- sprintf("python3 -c \"import scanpy as sc; ad=sc.read_h5ad('%s'); y=ad[:,'NDUFB7'].X.toarray().flatten() if hasattr(ad[:,'NDUFB7'].X,'toarray') else ad[:,'NDUFB7'].X.flatten(); import pandas as pd; pd.DataFrame({'NDUFB7':y}).to_csv('%s', index=False)\" 2>/dev/null", h5_candidates[1], tmpcsv)
    system(cmd)
    if (file.exists(tmpcsv)) {
      d <- fread(tmpcsv)
      if ("NDUFB7" %in% names(d)) {
        e121 <- as.numeric(d$NDUFB7); e121 <- e121[is.finite(e121)]
        if (length(e121) > 100) message("[PASS] h5ad: ", h5_candidates[1])
      }
    }
  }
}

# 路径3: RDS/Seurat
if (is.null(e121)) {
  rds_candidates <- list.files(c("01_data","03_results"), pattern = "GSE121893.*\\.rds", full.names = TRUE, recursive = TRUE)
  for (f in rds_candidates) {
    tryCatch({
      obj <- readRDS(f)
      if (inherits(obj, "Seurat") && "NDUFB7" %in% rownames(obj)) {
        e121 <- as.numeric(GetAssayData(obj, slot = "data")["NDUFB7", ])
        e121 <- e121[is.finite(e121)]
        if (length(e121) > 100) { message("[PASS] RDS: ", f); break }
      }
    }, error = function(e) NULL)
  }
}

if (is.null(e121) || length(e121) < 100) {
  message("[FAIL] GSE121893 NDUFB7 not available (need >100 cells, found ", ifelse(is.null(e121), 0, length(e121)), ")")
  cat("GSE121893_NOT_AVAILABLE\n", file = file.path(outdir, "V152_status.txt"))
  quit(save="no", status=1)
}

message("  n = ", length(e121), " | Zero% = ", round(mean(e121==0)*100, 1))

# === 2. BIC验证 ===
message("\n[2/3] BIC validation...")
x_nz <- e121[e121 > 0]
if (length(x_nz) < 50) {
  message("[FAIL] Too few non-zero values")
  quit(save="no", status=1)
}

bic_res <- data.frame()
for (G in 1:3) {
  tryCatch({
    mod <- densityMclust(x_nz, G = G, modelNames = "V", verbose = FALSE)
    bic_res <- rbind(bic_res, data.frame(Cohort = "GSE121893", G = G, BIC = mod$bic, logL = mod$loglik, n_nz = length(x_nz), stringsAsFactors = FALSE))
  }, error = function(e) message("  G=", G, " failed: ", conditionMessage(e)))
}

if (nrow(bic_res) > 0) {
  bic_res$Best <- bic_res$BIC == min(bic_res$BIC)
  fwrite(bic_res, file.path(outdir, "V152_GSE121893_BIC.csv"))
  message("\n=== GSE121893 BIC ===")
  print(bic_res)
  
  best_G <- bic_res$G[bic_res$Best]
  message("\n[RESULT] GSE121893 BIC favors G = ", best_G)
  
  # 与GSE183852比较
  v140a <- "03_results/V140_BIC_Resolution/V140A_ZI_GMM_183852.csv"
  if (file.exists(v140a)) {
    z183 <- fread(v140a)
    message("\n--- Comparison with GSE183852 ---")
    message("  GSE183852 ZI-GMM best G: ", z183$G[z183$Best == TRUE])
    message("  GSE121893 best G: ", best_G)
    if (best_G >= 2) {
      message("  [PASS] Cross-cohort multi-modality confirmed")
    } else {
      message("  [WARN] GSE121893 unimodal — platform/disease heterogeneity")
    }
  }
}

# === 3. Dip Test ===
message("\n[3/3] Dip test...")
has_diptest <- requireNamespace("diptest", quietly = TRUE)
if (!has_diptest) {
  tryCatch({ install.packages("diptest", repos = "https://cloud.r-project.org/", quiet = TRUE); has_diptest <- TRUE }, error = function(e) message("diptest install failed"))
}

if (has_diptest) {
  xj <- e121 + rnorm(length(e121), 0, 0.001)
  dt <- diptest::dip.test(xj, simulate.p.value = TRUE, B = 2000)
  dip_res <- data.frame(Cohort = "GSE121893", N = length(e121), Zero_Pct = round(mean(e121==0)*100,1), Dip = round(dt$statistic, 4), Dip_P = format(dt$p.value, digits=2, scientific = TRUE), Reject_Unimodal = dt$p.value < 0.05, stringsAsFactors = FALSE)
  fwrite(dip_res, file.path(outdir, "V152_GSE121893_Dip.csv"))
  message("  Dip = ", round(dt$statistic, 4), ", p = ", format(dt$p.value, digits=2))
  if (dt$p.value < 0.05) message("  [PASS] Rejects unimodality") else message("  [INFO] Unimodality not rejected")
}

message("[DONE] V152: ", outdir)
