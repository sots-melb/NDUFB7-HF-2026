#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V168_PanDeath_Solidity/V168C_ACSL4_GPX4")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V168C: GSE59867 ACSL4/GPX4比值与HF预后")
message("========================================")

# === 1. 读取临床数据 ===
clin_file <- "01_data/02_clinical/GSE59867_clinical_final.csv"
if(!file.exists(clin_file)){
  clin_file <- "01_data/02_clinical/GSE59867_clinical_clean.csv"
}

if(!file.exists(clin_file)){
  message("[FAIL] Clinical data not found")
  quit(save="no", status=1)
}

clin <- fread(clin_file)
message("\n[1/4] Clinical data: ", nrow(clin), " rows")
message("  Columns: ", paste(names(clin), collapse = ", "))

# 找关键列
ntprobnp_col <- grep("NT.proBNP|NTproBNP|BNP|proBNP", names(clin), value = TRUE, ignore.case = TRUE)[1]
outcome_col <- grep("outcome|death|transplant|LVAD|event|composite", names(clin), value = TRUE, ignore.case = TRUE)[1]
group_col <- grep("group|condition|disease|hf|status", names(clin), value = TRUE, ignore.case = TRUE)[1]

message("  NT-proBNP column: ", ifelse(is.na(ntprobnp_col), "NOT FOUND", ntprobnp_col))
message("  Outcome column: ", ifelse(is.na(outcome_col), "NOT FOUND", outcome_col))
message("  Group column: ", ifelse(is.na(group_col), "NOT FOUND", group_col))

# === 2. 读取/提取GSE59867表达数据中的ACSL4和GPX4 ===
# 尝试多个来源
expr_sources <- c(
  "03_results/02_figures_tables/GSE59867_NDUFB7_final.csv",
  "03_results/02_figures_tables/GSE59867_NDUFB7_expression_v2.csv",
  "03_results/02_figures_tables/GSE59867_NDUFB7_corrected.csv"
)

expr_df <- NULL
for(f in expr_sources){
  if(file.exists(f)){
    message("\n[2/4] 尝试读取表达数据: ", basename(f))
    df_temp <- fread(f)
    message("  Columns: ", paste(head(names(df_temp), 10), collapse = ", "))
    if(any(grepl("ACSL4|GPX4|NDUFB7", names(df_temp), ignore.case = TRUE))){
      expr_df <- df_temp
      message("  [PASS] 包含目标基因")
      break
    }
  }
}

# 如果上述文件没有ACSL4/GPX4，尝试从series matrix解析
if(is.null(expr_df) || !any(grepl("ACSL4", names(expr_df), ignore.case = TRUE))){
  message("\n[INFO] 现有文件不含ACSL4/GPX4，尝试从series matrix提取...")
  sm_file <- "01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz"
  if(file.exists(sm_file)){
    # series matrix可能有表达值
    lines <- readLines(sm_file, n = 100, warn = FALSE)
    has_matrix <- any(grepl("series_matrix_table", lines, ignore.case = TRUE))
    if(has_matrix){
      message("  Series matrix contains expression table, parsing...")
      # 找到matrix起始行
      # 这需要更复杂的解析，标记为需要手动处理
      message("  [WARN] Series matrix parsing complex — saving for manual extraction")
    }
  }
  
  # 回退：使用GPL570探针注释手动映射（如果用户有GPL570.annot.gz）
  gpl_file <- "01_data/01_raw_geo/GPL570/GPL570.annot.gz"
  if(!file.exists(gpl_file)) gpl_file <- "01_data/01_raw_geo/GPL570.annot.gz"
  if(!file.exists(gpl_file)) gpl_file <- "Downloads/GPL570.annot.gz"
  
  if(file.exists(gpl_file)){
    message("  GPL570 annotation found: ", gpl_file)
    message("  [ACTION] 需要手动提取ACSL4/GPX4探针表达 — 提供诊断信息")
  }
}

# === 3. 如果找到了ACSL4/GPX4，计算比值并分析 ===
if(!is.null(expr_df) && any(grepl("ACSL4", names(expr_df), ignore.case = TRUE)) && any(grepl("GPX4", names(expr_df), ignore.case = TRUE))){
  acsl4_col <- grep("ACSL4", names(expr_df), value = TRUE, ignore.case = TRUE)[1]
  gpx4_col <- grep("GPX4", names(expr_df), value = TRUE, ignore.case = TRUE)[1]
  
  expr_df$ACSL4_GPX4_Ratio <- as.numeric(expr_df[[acsl4_col]]) / (as.numeric(expr_df[[gpx4_col]]) + 0.001)
  expr_df$ACSL4_GPX4_Ratio <- log2(expr_df$ACSL4_GPX4_Ratio + 1)
  
  message("\n[3/4] ACSL4/GPX4比值统计:")
  message("  Mean: ", round(mean(expr_df$ACSL4_GPX4_Ratio, na.rm = TRUE), 3))
  message("  Range: ", round(min(expr_df$ACSL4_GPX4_Ratio, na.rm = TRUE), 3), " - ", round(max(expr_df$ACSL4_GPX4_Ratio, na.rm = TRUE), 3))
  
  # 合并临床数据
  merge_col <- intersect(names(expr_df), names(clin))[1]
  if(!is.na(merge_col)){
    merged <- merge(expr_df, clin, by = merge_col)
    message("  合并后样本: ", nrow(merged))
    
    # 与NT-proBNP相关
    if(!is.na(ntprobnp_col) && ntprobnp_col %in% names(merged)){
      np_vals <- as.numeric(merged[[ntprobnp_col]])
      ratio_vals <- merged$ACSL4_GPX4_Ratio
      valid <- is.finite(np_vals) & is.finite(ratio_vals)
      if(sum(valid) > 10){
        r <- cor(ratio_vals[valid], np_vals[valid], method = "spearman")
        message("\n  ACSL4/GPX4 vs NT-proBNP: rho = ", round(r, 3), " (n = ", sum(valid), ")")
        
        # 可视化
        df_plot <- data.frame(Ratio = ratio_vals[valid], NTproBNP = log2(np_vals[valid] + 1))
        p <- ggplot(df_plot, aes(x = Ratio, y = NTproBNP)) +
          geom_point(alpha = 0.5) + geom_smooth(method = "lm", se = TRUE) +
          labs(title = "ACSL4/GPX4 Ratio vs NT-proBNP (GSE59867)",
               subtitle = paste0("Spearman rho = ", round(r, 3)),
               x = "log2(ACSL4/GPX4 + 1)", y = "log2(NT-proBNP + 1)") +
          theme_minimal()
        ggsave(file.path(outdir, "V168C_ACSL4_GPX4_vs_NTproBNP.png"), p, width = 6, height = 5, dpi = 300)
      }
    }
    
    fwrite(merged, file.path(outdir, "V168C_ACSL4_GPX4_clinical_merged.csv"))
  }
} else {
  message("\n[WARN] ACSL4/GPX4 not found in available files")
  message("  [ACTION] 需要从GSE59867表达矩阵手动提取这两个基因")
  message("  建议: 使用GPL570探针注释映射，或从RAW.tar中重新处理CEL文件")
}

# 保存诊断信息
diag <- data.frame(
  Item = c("Clinical_File", "N_Clinical", "NTproBNP_Column", "Outcome_Column", "Expression_File", "ACSL4_Found", "GPX4_Found"),
  Value = c(basename(clin_file), nrow(clin), ifelse(is.na(ntprobnp_col), "NO", ntprobnp_col),
            ifelse(is.na(outcome_col), "NO", outcome_col), 
            ifelse(is.null(expr_df), "NO", "YES"),
            ifelse(!is.null(expr_df) && any(grepl("ACSL4", names(expr_df))), "YES", "NO"),
            ifelse(!is.null(expr_df) && any(grepl("GPX4", names(expr_df))), "YES", "NO")),
  stringsAsFactors = FALSE
)
fwrite(diag, file.path(outdir, "V168C_diagnosis.csv"))

message("\n[DONE] V168C: ", outdir)
