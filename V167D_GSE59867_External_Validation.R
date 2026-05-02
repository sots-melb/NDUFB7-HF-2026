#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
outdir <- file.path(PROJECT, "03_results/V167_PanDeath_Validation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V167D: GSE59867外部验证")
message("========================================")

# GSE59867表达数据
expr_file <- "03_results/02_figures_tables/GSE59867_NDUFB7_final.csv"
if(file.exists(expr_file)){
  df <- fread(expr_file)
  message("GSE59867 data: ", nrow(df), " rows")
  message("Columns: ", paste(names(df), collapse = ", "))
  
  # 检查是否有通路评分
  death_cols <- grep("Apoptosis|Necroptosis|Ferroptosis|Autophagy|Pyroptosis", names(df), value = TRUE)
  if(length(death_cols) > 0){
    message("\nDeath pathway columns found: ", paste(death_cols, collapse = ", "))
    
    ndufb7_col <- grep("NDUFB7", names(df), value = TRUE)[1]
    if(!is.na(ndufb7_col)){
      results <- data.frame()
      for(col in death_cols){
        r <- cor(df[[ndufb7_col]], df[[col]], method = "spearman", use = "pairwise.complete.obs")
        results <- rbind(results, data.frame(
          Cohort = "GSE59867",
          Pathway = col,
          Spearman_Rho = round(r, 3),
          N = sum(is.finite(df[[ndufb7_col]]) & is.finite(df[[col]])),
          stringsAsFactors = FALSE
        ))
      }
      
      fwrite(results, file.path(outdir, "V167D_GSE59867_external_validation.csv"))
      message("\n=== GSE59867外部验证 ===")
      print(results[order(abs(Spearman_Rho), decreasing = TRUE)])
      
      # 判定：是否与GSE57338方向一致？
      v151 <- fread("03_results/V151_C2_Discriminant/V151_discriminant_validation.csv")
      merged <- merge(results, v151, by.x = "Pathway", by.y = "Pathway", all.x = TRUE)
      merged$Direction_Consistent <- sign(merged$Spearman_Rho.x) == sign(merged$Spearman_Rho.y)
      message("\n方向一致性:")
      print(merged[, .(Pathway, Spearman_Rho.x, Spearman_Rho.y, Direction_Consistent)])
      
      n_consistent <- sum(merged$Direction_Consistent, na.rm = TRUE)
      message("\n[", ifelse(n_consistent >= 4, "PASS", "WARN"), "] ", n_consistent, "/", nrow(merged), " pathways directionally consistent across cohorts")
    }
  } else {
    message("[INFO] No death pathway scores in GSE59867 — need to compute from expression matrix")
    message("  标记为Revision任务")
  }
} else {
  message("[FAIL] GSE59867 file not found")
}

message("\n[DONE] V167D: ", outdir)
