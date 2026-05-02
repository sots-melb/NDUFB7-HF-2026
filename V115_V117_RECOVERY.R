#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mclust)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V115_V117_Recovery")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V115_V117: 跨模型数据修复与双峰检验")
message("========================================")

# 统一返回格式（必须与V114完全一致，否则V119 rbind失败）
test_bimodal <- function(expr_vec, label) {
  expr_vec <- as.numeric(expr_vec); expr_vec <- expr_vec[!is.na(expr_vec)]
  n <- length(expr_vec)
  res <- data.frame(label=label, n=n, g=NA, peak1=NA, peak2=NA, gap=NA, sd=NA, is_bimodal=NA, zero_pct=NA, mid_pct=NA, all_or_none=NA, stringsAsFactors=FALSE)
  if (n < 20) return(res)
  mod <- tryCatch(densityMclust(expr_vec, G=1:3, verbose=FALSE), error=function(e) NULL)
  if (is.null(mod)) {
    res$sd <- sd(expr_vec); res$zero_pct <- mean(expr_vec==0)*100; return(res)
  }
  best_g <- mod$G; peaks <- if(best_g>=2) sort(mod$parameters$mean) else NA
  gap <- if(best_g>=2 && length(peaks)>=2) abs(peaks[2]-peaks[1]) else 0
  sd_expr <- sd(expr_vec)
  is_bimodal <- best_g==2 && gap > 0.5*sd_expr
  zero_pct <- mean(expr_vec==0)*100
  nonzero <- expr_vec[expr_vec>0]
  high_thresh <- ifelse(length(nonzero)>0, quantile(nonzero,0.9), 0)
  mid_pct <- mean(expr_vec>0 & expr_vec<high_thresh)*100
  all_or_none <- zero_pct + mean(expr_vec>=high_thresh)*100
  data.frame(label=label, n=n, g=best_g, peak1=ifelse(length(peaks)>=1,peaks[1],NA), peak2=ifelse(length(peaks)>=2,peaks[2],NA), gap=gap, sd=sd_expr, is_bimodal=is_bimodal, zero_pct=zero_pct, mid_pct=mid_pct, all_or_none=all_or_none, stringsAsFactors=FALSE)
}

res_list <- list()

# --- A. GSE55296: 修复基因名大小写/空格问题 ---
f55296 <- "~/Downloads/GSE55296_count_data.txt.gz"
if (file.exists(f55296)) {
  message("[LOAD] GSE55296 count data...")
  dt <- fread(f55296, header=TRUE)
  gene_col <- colnames(dt)[1]
  genes <- trimws(toupper(as.character(dt[[gene_col]])))
  idx <- which(genes == "NDUFB7")[1]
  if (is.na(idx)) idx <- which(grepl("^NDUFB7$", genes))[1]
  if (!is.na(idx)) {
    expr <- suppressWarnings(as.numeric(dt[idx, -1, with=FALSE]))
    expr <- expr[!is.na(expr)]
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE55296_Bulk")
    message("[PASS] GSE55296 NDUFB7: n=", length(expr), " mean=", round(mean(expr),2))
  } else { message("[WARN] NDUFB7 not found in GSE55296 (checked ", length(genes), " rows)") }
} else { message("[SKIP] GSE55296 not found") }

# --- B. GSE121893 ---
f121893 <- "~/Downloads/GSE121893_human_heart_sc_umi.csv.gz"
if (file.exists(f121893)) {
  message("[LOAD] GSE121893...")
  dt <- fread(f121893, header=TRUE, nrows=50000)
  gene_col <- colnames(dt)[1]
  idx <- which(toupper(dt[[gene_col]])=="NDUFB7")[1]
  if (!is.na(idx)) {
    expr <- as.numeric(dt[idx, -1, with=FALSE])
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE121893_sc")
    message("[PASS] GSE121893: n=", length(expr))
  }
} else { message("[SKIP] GSE121893 not found") }

# --- C. GSE275031: 修复gzipped rds ---
f275031 <- "~/Downloads/GSE275031_integrated_seurat_obj.rds.gz"
if (file.exists(f275031)) {
  message("[LOAD] GSE275031 gzipped rds...")
  obj <- tryCatch(readRDS(gzfile(f275031)), error=function(e) NULL)
  if (!is.null(obj) && inherits(obj, "Seurat")) {
    gene <- ifelse("Ndufb7" %in% rownames(obj), "Ndufb7", ifelse("NDUFB7" %in% rownames(obj), "NDUFB7", NA))
    if (!is.na(gene)) {
      expr <- as.numeric(FetchData(obj, vars=gene)[[1]])
      res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE275031_Mouse_HFpEF")
      message("[PASS] GSE275031: n=", length(expr))
    } else { message("[WARN] Ndufb7/NDUFB7 not in GSE275031 rownames") }
  } else { message("[WARN] GSE275031 not readable as Seurat (may need gunzip first)") }
} else { message("[SKIP] GSE275031 not found") }

# --- D. GSE168742: 修复非Seurat对象 ---
f168 <- list.files("01_data", pattern="GSE168742.*\\.rds$|168742.*human.*\\.rds$", recursive=TRUE, full.names=TRUE)[1]
if (is.na(f168)) f168 <- list.files("03_results", pattern="GSE168742.*\\.rds$", recursive=TRUE, full.names=TRUE)[1]
if (!is.na(f168) && file.exists(f168)) {
  message("[LOAD] GSE168742: ", basename(f168))
  obj <- readRDS(f168)
  if (inherits(obj, "Seurat")) {
    expr <- as.numeric(FetchData(obj, vars="NDUFB7")$NDUFB7)
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE168742_sc")
    message("[PASS] GSE168742 Seurat: n=", length(expr))
  } else if (is.matrix(obj) || is.data.frame(obj)) {
    rnames <- rownames(obj)
    if (!is.null(rnames) && "NDUFB7" %in% rnames) {
      expr <- as.numeric(obj["NDUFB7", ])
      res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE168742_matrix")
      message("[PASS] GSE168742 matrix (row=genes): n=", length(expr))
    } else {
      cnames <- colnames(obj)
      if (!is.null(cnames) && "NDUFB7" %in% cnames) {
        expr <- as.numeric(obj[, "NDUFB7"])
        res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE168742_matrix_col")
        message("[PASS] GSE168742 matrix (col=genes): n=", length(expr))
      } else { message("[WARN] NDUFB7 not in dimnames") }
    }
  } else { message("[WARN] Unknown object type: ", paste(class(obj), collapse="/")) }
} else { message("[SKIP] GSE168742 not found") }

# --- E. GSE106118 胚胎 ---
f106118 <- "~/Downloads/GSE106118_UMI_count_merge.txt.gz"
if (file.exists(f106118)) {
  message("[LOAD] GSE106118...")
  dt <- fread(f106118, header=TRUE, nrows=50000)
  gene_col <- colnames(dt)[1]
  idx <- which(toupper(dt[[gene_col]])=="NDUFB7")[1]
  if (!is.na(idx)) {
    expr <- as.numeric(dt[idx, -1, with=FALSE])
    res_list[[length(res_list)+1]] <- test_bimodal(expr, "GSE106118_Embryo")
    message("[PASS] GSE106118: n=", length(expr))
  }
} else { message("[SKIP] GSE106118 not found") }

# --- F. GSE116250 RPKM ---
f116250 <- list.files(c("~/Downloads","01_data"), pattern="GSE116250.*rpkm.*\\.txt\\.gz$|GSE116250.*\\.rds$", full.names=TRUE)[1]
if (!is.na(f116250) && file.exists(f116250)) {
  message("[LOAD] GSE116250: ", basename(f116250))
  if (grepl("\\.rds$", f116250)) {
    obj <- readRDS(f116250)
    # 灵活处理
    message("[INFO] GSE116250 rds loaded, type=", paste(class(obj), collapse="/"))
  } else {
    dt <- fread(f116250, header=TRUE, nrows=100)
    message("[INFO] GSE116250 txt loaded, cols=", ncol(dt))
  }
} else { message("[SKIP] GSE116250 not found") }

# --- G. GSE141910 ---
f141910 <- list.files(c("~/Downloads","01_data"), pattern="GSE141910", full.names=TRUE)[1]
if (!is.na(f141910) && file.exists(f141910)) {
  message("[INFO] GSE141910 exists: ", basename(f141910), " (format-specific extraction needed)")
} else { message("[SKIP] GSE141910 not found") }

# --- 保存 ---
if (length(res_list)>0) {
  df_res <- do.call(rbind, res_list)
  write.csv(df_res, file.path(outdir,"V115_V117_recovered_bimodal.csv"), row.names=FALSE)
  message("\n=== 修复后跨模型双峰汇总 ===")
  print(df_res[,c("label","n","g","is_bimodal","zero_pct","mid_pct","all_or_none")])
  
  p <- ggplot(df_res, aes(x=reorder(label, all_or_none), y=all_or_none, fill=is_bimodal)) +
    geom_bar(stat="identity") + coord_flip() +
    scale_fill_manual(values=c("TRUE"="#FDE725","FALSE"="#440154"), na.value="grey50") +
    labs(title="All-or-None Index (Recovered Data)", x="", y="Zero% + Top10%") + theme_minimal()
  ggsave(file.path(outdir,"V115_cross_model_bar.png"), p, width=8, height=5, dpi=300)
} else { message("[WARN] No data recovered") }

message("[DONE] V115_V117: ", outdir)
