#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)
library(metafor)

cat("\n========== Phase 2 Pillar 1: NDUFB7荟萃分析 ==========\n")

# 数据集配置
datasets <- list(
  list(id = "GSE57338",  file = "01_data/01_raw_geo/GSE57338/GSE57338_gene_level.rds", 
       ndufb7 = "NDUFB7", platform = "Affy ST"),
  list(id = "GSE141910", file = "01_data/01_raw_geo/GSE141910/GSE141910_merged_matrix.rds", 
       ndufb7 = "ENSG00000167996", platform = "RNA-seq"),
  list(id = "GSE5406",   file = "01_data/01_raw_geo/GSE5406/GSE5406_gene_level.rds", 
       ndufb7 = "NDUFB7", platform = "U133A"),
  list(id = "GSE79962",  file = "01_data/01_raw_geo/GSE79962/GSE79962_gene_level.rds", 
       ndufb7 = "NDUFB7", platform = "Affy ST"),
  list(id = "GSE116250", file = "01_data/01_raw_geo/GSE116250/GSE116250_rpkm_matrix.rds", 
       ndufb7 = "ENSG00000167996", platform = "RNA-seq (RPKM)")
)

# 分组提取（简化版，基于已知模式）
get_groups <- function(gse_id, samples) {
  n <- length(samples)
  if (gse_id == "GSE116250") {
    return(ifelse(grepl("^NF", samples), "Control", "HF"))
  }
  # 其他数据集从Series Matrix提取（简化处理）
  sm_file <- paste0("01_data/01_raw_geo/", gse_id, "/", gse_id, "_series_matrix.txt.gz")
  if (file.exists(sm_file)) {
    con <- gzfile(sm_file, "r")
    lines <- readLines(con)
    close(con)
    meta <- lines[1:min(100, length(lines))]
    for (field in c("source_name", "title", "characteristics_ch1")) {
      fld <- meta[grepl(paste0("!Sample_", field), meta, ignore.case = TRUE)]
      if (length(fld) > 0) {
        vals <- gsub('^"|"$', '', unlist(strsplit(fld[1], "\t"))[-1])
        vals_l <- tolower(vals)
        grp <- rep(NA, n)
        grp[grepl("failing|dcm|icm|hcm|cardiomyopathy|patient|disease|hf$", vals_l)] <- "HF"
        grp[grepl("non.failing|nonfailing|control|normal|healthy|donor|nf$", vals_l)] <- "Control"
        if (sum(!is.na(grp)) > n * 0.3) return(grp)
      }
    }
  }
  return(rep(NA, n))
}

# Cohen's d计算
cohens_d <- function(x, y) {
  nx <- sum(!is.na(x)); ny <- sum(!is.na(y))
  if (nx < 2 || ny < 2) return(c(d = NA, v = NA))
  mx <- mean(x, na.rm = TRUE); my <- mean(y, na.rm = TRUE)
  sp <- sqrt(((nx-1)*sd(x,na.rm=T)^2 + (ny-1)*sd(y,na.rm=T)^2) / (nx+ny-2))
  d <- (mx - my) / sp
  v <- (nx+ny)/(nx*ny) + d^2/(2*(nx+ny))
  c(d = d, v = v)
}

# 处理每个数据集
results <- list()
for (ds in datasets) {
  cat("\n[处理]", ds$id, "...\n")
  if (!file.exists(ds$file)) { cat("  ❌ 文件缺失\n"); next }
  
  obj <- readRDS(ds$file)
  mat <- if ("exprs" %in% names(obj)) obj$exprs else obj
  
  if (!ds$ndufb7 %in% rownames(mat)) { cat("  ❌ NDUFB7缺失\n"); next }
  
  expr <- as.numeric(mat[ds$ndufb7, ])
  samples <- colnames(mat)
  groups <- get_groups(ds$id, samples)
  
  hf <- expr[groups == "HF"]
  ctrl <- expr[groups == "Control"]
  
  if (length(hf) < 3 || length(ctrl) < 3) {
    cat("  ⚠️ 分组不足 (HF=", length(hf), ", Ctrl=", length(ctrl), ")\n")
    next
  }
  
  # Z-score标准化
  expr_z <- as.vector(scale(expr))
  eff <- cohens_d(expr_z[groups == "HF"], expr_z[groups == "Control"])
  
  cat("  ✅ d =", round(eff["d"], 3), "| HF=", length(hf), "Ctrl=", length(ctrl), "\n")
  
  results[[ds$id]] <- c(yi = eff["d"], vi = eff["v"], n = length(expr))
}

# 荟萃分析
if (length(results) >= 2) {
  meta_df <- do.call(rbind, lapply(names(results), function(x) {
    data.frame(id = x, yi = results[[x]]["yi"], vi = results[[x]]["vi"], n = results[[x]]["n"])
  }))
  
  cat("\n效应量汇总:\n")
  print(meta_df)
  
  res <- rma(yi = yi, vi = vi, data = meta_df, method = "REML",
             slab = paste0(id, " (n=", n, ")"))
  
  cat("\n========== 随机效应模型结果 ==========\n")
  print(res)
  
  # 森林图
  dir.create("03_results/07_figures", showWarnings = FALSE, recursive = TRUE)
  pdf("03_results/07_figures/Fig1B_forest.pdf", width = 11, height = 7)
  forest(res, main = "NDUFB7 in Heart Failure - Meta-Analysis",
         xlab = "Standardized Mean Difference (Cohen's d)",
         mlab = paste0("RE Model (I² = ", round(res$I2, 1), "%)"))
  dev.off()
  
  png("03_results/07_figures/Fig1B_forest.png", width = 1000, height = 700, res = 150)
  forest(res, main = "NDUFB7 in Heart Failure",
         xlab = "SMD (Cohen's d)",
         mlab = paste0("RE Model (I² = ", round(res$I2, 1), "%)"))
  dev.off()
  
  cat("\n[保存] Figure 1B 已生成\n")
  cat("  PDF: 03_results/07_figures/Fig1B_forest.pdf\n")
  cat("  PNG: 03_results/07_figures/Fig1B_forest.png\n")
  cat("\n汇总效应: d =", round(res$b, 3), "[", round(res$ci.lb, 3), ",", round(res$ci.ub, 3), "]\n")
  cat("p值:", format(res$pval, digits = 3), "| I²:", round(res$I2, 1), "%\n")
  
  # Egger检验
  if (nrow(meta_df) >= 3) {
    eg <- regtest(res)
    cat("Egger检验 p =", round(eg$pval, 4), "\n")
  }
} else {
  cat("❌ 有效数据集不足\n")
}
