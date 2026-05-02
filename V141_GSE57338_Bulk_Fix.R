#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(ggplot2); library(GEOquery); library(data.table) })
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V141_GSE57338_Bulk_Fix")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V141: GSE57338 Bulk铁死亡评分（修复版）")
message("========================================")

# ==================== 多路径加载GSE57338 ====================
exprs <- NULL
sample_names <- NULL

# 路径1: 搜索任何GSE57338 RDS
rds_files <- list.files(c("01_data","03_results"), pattern = "GSE57338.*\\.rds$", full.names = TRUE, recursive = TRUE)
if (length(rds_files) > 0) {
  for (f in rds_files) {
    tryCatch({
      obj <- readRDS(f)
      if (is.list(obj) && "exprs" %in% names(obj)) {
        exprs <- obj$exprs
        sample_names <- colnames(exprs)
        message("[RDS] Loaded from ", basename(f), ": ", nrow(exprs), "x", ncol(exprs))
        break
      } else if (is.matrix(obj) || is.data.frame(obj)) {
        exprs <- as.matrix(obj)
        sample_names <- colnames(exprs)
        message("[RDS] Loaded matrix from ", basename(f), ": ", nrow(exprs), "x", ncol(exprs))
        break
      }
    }, error = function(e) NULL)
  }
}

# 路径2: 搜索V135 CSV输出
if (is.null(exprs)) {
  csv_files <- list.files("03_results", pattern = "V135.*GSE57338.*\\.csv", full.names = TRUE, recursive = TRUE)
  if (length(csv_files) > 0) {
    for (f in csv_files) {
      tryCatch({
        d <- fread(f)
        if (all(c("Sample","NDUFB7") %in% names(d))) {
          # 这是样本级别的汇总，不是矩阵，无法做通路评分
          message("[CSV] Found V135 summary but not expression matrix: ", basename(f))
        }
      }, error = function(e) NULL)
    }
  }
}

# 路径3: 在线GEO加载（最终回退）
if (is.null(exprs)) {
  message("[GEO] Attempting online download of GSE57338...")
  tryCatch({
    gse <- getGEO("GSE57338", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
    exprs <- exprs(gse)
    sample_names <- colnames(exprs)
    message("[GEO] Loaded from GEO online: ", nrow(exprs), "x", ncol(exprs))
  }, error = function(e) {
    message("[FAIL] Online GEO failed: ", conditionMessage(e))
  })
}

if (is.null(exprs)) stop("GSE57338 expression matrix unavailable. Cannot proceed.")

# 确认NDUFB7存在
if (!("NDUFB7" %in% rownames(exprs))) {
  # 尝试探针映射
  message("[WARN] NDUFB7 not in rownames. Attempting probe-to-gene mapping...")
  # 简化：如果行名是探针ID，尝试从GPL11532映射
  stop("NDUFB7 not found in expression matrix rownames. Manual probe mapping required.")
}

ndufb7_expr <- as.numeric(exprs["NDUFB7", ])
names(ndufb7_expr) <- sample_names
message("[PASS] NDUFB7 mean: ", round(mean(ndufb7_expr, na.rm=TRUE), 2))

# ==================== 铁死亡评分 ====================
ferro_genes <- c("FTL","SAT1","NFE2L2","SLC7A11","ACSL4","GPX4","FTH1","TFRC","PTGS2","KEAP1")
avail <- intersect(ferro_genes, rownames(exprs))
message("[1/4] Ferroptosis genes available: ", length(avail), "/", length(ferro_genes), " | ", paste(avail, collapse=","))

if (length(avail) >= 3) {
  ferro_mat <- exprs[avail, , drop = FALSE]
  
  # 安全z-score（逐基因）
  ferro_z_list <- lapply(1:nrow(ferro_mat), function(i) {
    g <- as.numeric(ferro_mat[i, ])
    m <- mean(g, na.rm = TRUE); s <- sd(g, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(g)))
    (g - m) / s
  })
  ferro_z <- do.call(rbind, ferro_z_list)
  colnames(ferro_z) <- sample_names
  
  ferro_score <- colMeans(ferro_z, na.rm = TRUE)
  names(ferro_score) <- sample_names
  
  # 对齐样本
  common <- intersect(names(ndufb7_expr), names(ferro_score))
  if (length(common) == 0) {
    message("[WARN] Name mismatch, trying positional alignment")
    if (length(ndufb7_expr) == length(ferro_score)) {
      common <- seq_along(ndufb7_expr)
      names(ferro_score) <- names(ndufb7_expr)
      common <- intersect(names(ndufb7_expr), names(ferro_score))
    }
  }
  
  ndufb7_al <- ndufb7_expr[common]
  ferro_al <- ferro_score[common]
  
  valid <- is.finite(ndufb7_al) & is.finite(ferro_al)
  ndufb7_clean <- ndufb7_al[valid]
  ferro_clean <- ferro_al[valid]
  
  message("  Aligned samples: ", length(ndufb7_clean))
  
  # Spearman
  ct <- cor.test(ndufb7_clean, ferro_clean, method = "spearman")
  message("  NDUFB7 vs Ferroptosis: rho=", round(ct$estimate, 3), 
          " p=", format(ct$p.value, digits=2, scientific=TRUE))
  
  # 分组比较
  med <- median(ndufb7_clean)
  grp <- ifelse(ndufb7_clean < med, "NDUFB7_Low", "NDUFB7_High")
  wt <- wilcox.test(ferro_clean ~ grp)
  message("  Low vs High Wilcoxon: p=", format(wt$p.value, digits=2, scientific=TRUE))
  
  # 保存
  df <- data.frame(
    sample = names(ndufb7_clean),
    NDUFB7 = ndufb7_clean,
    ferroptosis_score = ferro_clean,
    NDUFB7_group = grp,
    stringsAsFactors = FALSE
  )
  write.csv(df, file.path(outdir, "V141_bulk_ferroptosis.csv"), row.names = FALSE)
  
  # 可视化
  p <- ggplot(df, aes(x = NDUFB7_group, y = ferroptosis_score, fill = NDUFB7_group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("NDUFB7_Low"="#440154", "NDUFB7_High"="#35B779")) +
    labs(title = "Bulk Ferroptosis Score by NDUFB7 (GSE57338)",
         subtitle = paste0("Spearman rho=", round(ct$estimate, 3), " p=", signif(ct$p.value, 2),
                          " | Wilcoxon p=", signif(wt$p.value, 2)),
         x = "NDUFB7 Group", y = "Ferroptosis Score (z-mean)") +
    theme_minimal(base_size = 12)
  ggsave(file.path(outdir, "V141_bulk_ferroptosis_boxplot.png"), p, width = 6, height = 5, dpi = 300)
  
  # ACSL4/GPX4比值
  if (all(c("ACSL4","GPX4") %in% rownames(exprs))) {
    ratio <- as.numeric(exprs["ACSL4", ]) / (as.numeric(exprs["GPX4", ]) + 1e-6)
    names(ratio) <- sample_names
    ratio_al <- ratio[names(ndufb7_clean)]
    valid2 <- is.finite(ratio_al)
    if (sum(valid2) > 10) {
      cr <- cor.test(ndufb7_clean[valid2], ratio_al[valid2], method = "spearman")
      message("  ACSL4/GPX4 vs NDUFB7: rho=", round(cr$estimate, 3),
              " p=", format(cr$p.value, digits=2, scientific=TRUE))
      df$ACSL4_GPX4_ratio <- ratio_al
      write.csv(df, file.path(outdir, "V141_bulk_ferroptosis.csv"), row.names = FALSE)
    }
  }
  
  message("[DONE] V141: ", outdir)
} else {
  message("[FAIL] Insufficient ferroptosis genes available (need >=3, found ", length(avail), ")")
}
