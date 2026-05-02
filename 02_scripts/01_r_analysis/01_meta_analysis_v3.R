#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)
library(metafor)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     Phase 2 Pillar 1: NDUFB7荟萃分析 v3（最终修复版）    ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# ==================== 1. 数据集配置 ====================
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

# ==================== 2. 终极分组提取函数 ====================
extract_groups_final <- function(gse_id, sample_names) {
  n <- length(sample_names)
  groups <- rep(NA_character_, n)
  
  # GSE116250: NF=Control, DCM/ICM=HF
  if (gse_id == "GSE116250") {
    groups <- ifelse(grepl("^NF", sample_names), "Control", "HF")
    cat("  [分组] GSE116250: NF=Control, DCM/ICM=HF\n")
    return(groups)
  }
  
  # GSE141910: 从列名推断（基于排查结果：包含DCM/NF/HCM等）
  if (gse_id == "GSE141910") {
    s_lower <- tolower(sample_names)
    # 关键词匹配
    is_ctrl <- grepl("nf|non.?fail|normal|control|healthy|donor", s_lower)
    is_hf <- grepl("dcm|icm|hcm|fail|cardiomyopathy|patient|disease", s_lower)
    
    if (sum(is_ctrl | is_hf) > n * 0.5) {
      groups[is_ctrl] <- "Control"
      groups[is_hf] <- "HF"
      cat("  [分组] GSE141910: 从列名关键词推断\n")
    } else {
      # 备用：基于GEO描述（180 DCM + 186 NF）
      cat("  [分组] GSE141910: 使用GEO已知信息（180 HF + 186 Control）\n")
      groups[1:180] <- "HF"
      groups[181:366] <- "Control"
    }
    return(groups)
  }
  
  # 其他数据集：从Series Matrix提取
  sm_file <- paste0("01_data/01_raw_geo/", gse_id, "/", gse_id, "_series_matrix.txt.gz")
  if (file.exists(sm_file)) {
    con <- gzfile(sm_file, "r")
    lines <- readLines(con)
    close(con)
    
    begin <- grep("^!series_matrix_table_begin", lines, ignore.case = TRUE)
    meta_lines <- if (length(begin) > 0) lines[1:(begin[1]-1)] else lines
    
    # 按优先级查找
    for (field in c("^!Sample_source_name", "^!Sample_title", "^!Sample_characteristics_ch1")) {
      fld_lines <- meta_lines[grepl(field, meta_lines, ignore.case = TRUE)]
      if (length(fld_lines) > 0) {
        vals <- gsub('^"|"$', '', unlist(strsplit(fld_lines[1], "\t"))[-1])
        vals_lower <- tolower(vals)
        
        is_hf <- grepl("failing|dcm|icm|hcm|dilated|ischemic|hypertrophic|cardiomyopathy|hf$|heart failure|patient|disease", vals_lower)
        is_ctrl <- grepl("non.failing|nonfailing|^nf$|control|normal|healthy|donor|^ctrl$", vals_lower)
        
        groups[is_hf] <- "HF"
        groups[is_ctrl] <- "Control"
        
        if (sum(!is.na(groups)) > n * 0.3) break
      }
    }
  }
  
  return(groups)
}

# ==================== 3. Cohen's d计算 ====================
cohens_d <- function(x, y) {
  nx <- sum(!is.na(x)); ny <- sum(!is.na(y))
  if (nx < 2 || ny < 2) return(c(d = NA, v = NA))
  mx <- mean(x, na.rm = TRUE); my <- mean(y, na.rm = TRUE)
  sp <- sqrt(((nx-1)*sd(x,na.rm=T)^2 + (ny-1)*sd(y,na.rm=T)^2) / (nx+ny-2))
  d <- (mx - my) / sp
  v <- (nx+ny)/(nx*ny) + d^2/(2*(nx+ny))
  c(d = d, v = v)
}

# ==================== 4. 逐数据集处理 ====================
results <- list()
group_summary <- list()

for (ds in datasets) {
  cat("\n==========", ds$id, "==========\n")
  
  if (!file.exists(ds$file)) { cat("  ❌ 文件不存在\n"); next }
  
  obj <- readRDS(ds$file)
  mat <- if ("exprs" %in% names(obj)) obj$exprs else obj
  
  if (!ds$ndufb7 %in% rownames(mat)) { cat("  ❌ NDUFB7缺失\n"); next }
  
  expr_raw <- as.numeric(mat[ds$ndufb7, ])
  sample_names <- colnames(mat)
  
  cat("  原始表达范围:", paste(round(range(expr_raw, na.rm = TRUE), 2), collapse = " ~ "), "\n")
  
  # 提取分组
  groups <- extract_groups_final(ds$id, sample_names)
  
  cat("  分组统计:\n")
  print(table(groups, useNA = "ifany"))
  
  # 保存分组
  group_summary[[ds$id]] <- data.frame(
    Dataset = ds$id, Sample = sample_names, Group = groups,
    NDUFB7_Raw = expr_raw, stringsAsFactors = FALSE
  )
  
  hf_idx <- which(groups == "HF")
  ctrl_idx <- which(groups == "Control")
  
  if (length(hf_idx) < 3 || length(ctrl_idx) < 3) {
    cat("  ⚠️ 分组样本不足，跳过\n")
    next
  }
  
  # 样本失衡警告
  ratio <- length(hf_idx) / length(ctrl_idx)
  if (ratio > 5 || ratio < 0.2) {
    cat("  ⚠️ 样本失衡（HF:Control =", length(hf_idx), ":", length(ctrl_idx), "）\n")
  }
  
  # Z-score标准化
  expr_z <- as.vector(scale(expr_raw))
  
  # 效应量（HF vs Control）
  eff <- cohens_d(expr_z[hf_idx], expr_z[ctrl_idx])
  
  cat("  ✅ Cohen's d =", round(eff["d"], 3), "| Var =", round(eff["v"], 4), "\n")
  cat("     HF n=", length(hf_idx), ", Control n=", length(ctrl_idx), "\n")
  
  results[[ds$id]] <- list(
    id = ds$id, platform = ds$platform,
    yi = eff["d"], vi = eff["v"],
    n_hf = length(hf_idx), n_ctrl = length(ctrl_idx),
    n_total = length(expr_raw)
  )
}

# 保存分组
dir.create("03_results/08_tables", showWarnings = FALSE, recursive = TRUE)
group_df <- do.call(rbind, group_summary)
write.csv(group_df, "03_results/08_tables/TableS1_sample_groups_v3.csv", row.names = FALSE)
cat("\n[保存] TableS1_sample_groups_v3.csv\n")

# ==================== 5. 荟萃分析 ====================
if (length(results) < 2) {
  cat("\n❌ 有效数据集不足2个\n")
  quit(status = 1)
}

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║              metafor 随机效应模型                          ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

meta_df <- do.call(rbind, lapply(names(results), function(x) {
  data.frame(id = x, yi = results[[x]]$yi, vi = results[[x]]$vi, 
             n = results[[x]]$n_total, stringsAsFactors = FALSE)
}))

cat("\n效应量汇总:\n")
print(meta_df)

res <- rma(yi = yi, vi = vi, data = meta_df, method = "REML",
           slab = paste0(id, " (n=", n, ")"))

cat("\n========== 荟萃分析结果 ==========\n")
print(res)

# ==================== 6. 森林图 ====================
dir.create("03_results/07_figures", showWarnings = FALSE, recursive = TRUE)

pdf("03_results/07_figures/Fig1B_forest_v3.pdf", width = 13, height = 9)
forest(res, main = "NDUFB7 Expression in Heart Failure: Meta-Analysis (k=5)",
       xlab = "Standardized Mean Difference (Cohen's d)",
       mlab = paste0("RE Model (I² = ", round(res$I2, 1), "%, p = ", format(res$QEp, digits = 3), ")"),
       cex = 0.9, col = "darkblue", border = "darkblue", addfit = TRUE, addpred = TRUE)
dev.off()

png("03_results/07_figures/Fig1B_forest_v3.png", width = 1200, height = 900, res = 150)
forest(res, main = "NDUFB7 in Heart Failure (k=5)",
       xlab = "SMD (Cohen's d)",
       mlab = paste0("RE Model (I² = ", round(res$I2, 1), "%)"))
dev.off()

cat("\n[保存] Figure 1B v3 已生成\n")

# ==================== 7. 诊断与解读 ====================
cat("\n========== 诊断统计量 ==========\n")
cat("纳入研究数: k =", length(results), "\n")
cat("异质性 I² =", round(res$I2, 1), "% (", 
    ifelse(res$I2 < 25, "低", ifelse(res$I2 < 50, "中等", "高")), "异质性)\n")
cat("Q统计量 p值 =", format(res$QEp, digits = 4), "\n")
cat("汇总效应量 d =", round(res$b[1], 3), "[", round(res$ci.lb, 3), ",", round(res$ci.ub, 3), "]\n")
cat("汇总效应 p值 =", format(res$pval, digits = 6), "\n")

# 各数据集方向一致性检查
cat("\n========== 方向一致性检查 ==========\n")
for (i in 1:nrow(meta_df)) {
  direction <- ifelse(meta_df$yi[i] > 0, "HF ↑", "HF ↓")
  cat(meta_df$id[i], ": d =", round(meta_df$yi[i], 3), direction, "\n")
}

# Egger检验
if (nrow(meta_df) >= 3) {
  cat("\n[Egger检验] 发表偏倚:\n")
  egger <- tryCatch(regtest(res, model = "lm"), error = function(e) NULL)
  if (!is.null(egger) && is.numeric(egger$est)) {
    cat("  p =", round(egger$pval, 4), "\n")
  } else {
    cat("  跳过（研究数不足或计算错误）\n")
  }
}

cat("\n========================================\n")
cat("🎉 Pillar 1 荟萃分析 v3 完成！\n")
cat("纳入数据集:", length(results), "/ 5\n")
cat("========================================\n")
