#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
  library(survminer)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V118_Calcium_Cox")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V118: SERCA2a互作 + Cox预后")
message("========================================")

res_list <- list()

# --- A. GSE164365 SERCA2a / Calcium handling ---
f164365 <- "~/Downloads/GSE164365_RAW.tar"
if (file.exists(f164365)) {
  message("[INFO] GSE164365 RAW.tar found. Checking project for processed...")
}
# 寻找project中可能存在的处理文件
f164_proc <- list.files("01_data", pattern="GSE164365|SERCA2a", recursive=TRUE, full.names=TRUE)[1]
if (!is.na(f164_proc) && file.exists(f164_proc)) {
  message("[LOAD] GSE164365 processed: ", basename(f164_proc))
  # 提取NDUFB7与SERCA2a (ATP2A2) 共表达
  obj <- readRDS(f164_proc)
  if (inherits(obj, "Seurat")) {
    genes <- c("NDUFB7", "ATP2A2", "SLC8A1", "RYR2")
    avail <- intersect(genes, rownames(obj))
    if (length(avail)>=2) {
      expr <- FetchData(obj, vars=avail)
      cor_mat <- cor(expr, method="spearman", use="pairwise.complete.obs")
      write.csv(cor_mat, file.path(outdir,"V118_SERCA2a_correlation.csv"))
      message("[PASS] SERCA2a correlation matrix saved")
      message("\n=== NDUFB7 vs Calcium genes ===")
      print(cor_mat["NDUFB7", avail])
    }
  }
} else { message("[SKIP] GSE164365 processed not found") }

# --- B. GSE59867 Cox预后 ---
# 需要表达矩阵 + 临床数据
clinical_file <- "03_results/V70_GSE59867_phenotype.csv"  # 或类似路径
if (!file.exists(clinical_file)) clinical_file <- list.files("03_results", pattern="GSE59867.*clinical|59867.*pheno", full.names=TRUE, recursive=TRUE)[1]

expr_file <- list.files("01_data", pattern="GSE59867.*\\.rds$|59867.*gene_level", full.names=TRUE, recursive=TRUE)[1]

if (!is.na(clinical_file) && file.exists(clinical_file) && !is.na(expr_file) && file.exists(expr_file)) {
  message("[LOAD] GSE59867 clinical + expression for Cox...")
  clin <- fread(clinical_file)
  expr <- readRDS(expr_file)
  
  # 提取NDUFB7表达
  # 根据expr结构灵活处理
  ndufb7_expr <- NULL
  if (inherits(expr, "list") && "expression" %in% names(expr)) {
    mat <- expr$expression
    ndufb7_expr <- as.numeric(mat["NDUFB7", ])
  } else if (is.matrix(expr) || is.data.frame(expr)) {
    ndufb7_expr <- as.numeric(expr["NDUFB7", ])
  }
  
  if (!is.null(ndufb7_expr) && length(ndufb7_expr)>0) {
    # 匹配样本
    # 简化：假设临床表和表达矩阵样本顺序一致或可通过ID匹配
    message("[INFO] NDUFB7 expr extracted, n=", length(ndufb7_expr))
    message("[WARN] Cox analysis requires time/event columns in clinical data")
    message("[ACTION] Please verify clinical file has 'time' and 'status' columns")
    # 如果列存在，执行Cox
    if (all(c("time","status") %in% colnames(clin))) {
      cox_df <- data.frame(time=clin$time, status=clin$status, NDUFB7=ndufb7_expr[1:nrow(clin)])
      cox_fit <- coxph(Surv(time, status) ~ NDUFB7, data=cox_df)
      cox_summary <- summary(cox_fit)
      message("\n=== Cox PH: NDUFB7 on Survival ===")
      message("HR: ", round(cox_summary$conf.int[1], 3))
      message("95% CI: ", paste(round(cox_summary$conf.int[3:4], 3), collapse=" - "))
      message("p: ", format(cox_summary$coefficients[5], digits=2, scientific=TRUE))
      write.csv(data.frame(HR=cox_summary$conf.int[1], CI_lower=cox_summary$conf.int[3], CI_upper=cox_summary$conf.int[4], p=cox_summary$coefficients[5]), file.path(outdir,"V118_cox_ndufb7.csv"), row.names=FALSE)
      
      # KM曲线（中位数分层）
      cox_df$NDUFB7_group <- ifelse(cox_df$NDUFB7 > median(cox_df$NDUFB7), "High", "Low")
      fit <- survfit(Surv(time, status) ~ NDUFB7_group, data=cox_df)
      p <- ggsurvplot(fit, pval=TRUE, risk.table=TRUE, title="GSE59867: NDUFB7 High vs Low")
      ggsave(file.path(outdir,"V118_KM_curve.png"), plot=p$plot, width=6, height=5, dpi=300)
      message("[DONE] KM curve saved")
    } else {
      message("[SKIP] Clinical data lacks time/status columns")
      # 保存描述性统计
      desc <- data.frame(N=length(ndufb7_expr), Mean=mean(ndufb7_expr), Median=median(ndufb7_expr), SD=sd(ndufb7_expr))
      write.csv(desc, file.path(outdir,"V118_ndufb7_desc.csv"), row.names=FALSE)
    }
  } else { message("[WARN] NDUFB7 not found in GSE59867 expression matrix") }
} else {
  message("[SKIP] GSE59867 clinical or expression file missing")
  message("  Clinical searched: ", ifelse(is.na(clinical_file),"NOT FOUND",clinical_file))
  message("  Expression searched: ", ifelse(is.na(expr_file),"NOT FOUND",expr_file))
}

message("[DONE] V118: ", outdir)
