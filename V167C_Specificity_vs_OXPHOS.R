#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
outdir <- file.path(PROJECT, "03_results/V167_PanDeath_Validation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V167C: NDUFB7特异性——vs其他线粒体基因")
message("========================================")

# 线粒体复合体I基因
complex1_genes <- c("NDUFB7", "NDUFA13", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11", "NDUFA1", "NDUFA2", "NDUFA3")

# 使用GSE57338基因水平数据
gene_file <- "03_results/01_seurat_objects/GSE57338_gene_level.rds"
if(file.exists(gene_file)){
  message("[INFO] GSE57338 gene level data found, but RDS读取需要特定环境")
  message("  回退：使用series matrix中的探针水平数据")
}

# 回退：从GSE57338 expression matrix提取
expr_file <- "03_results/02_figures_tables/GSE57338_ferroptosis_sample_scores_v3.csv"
if(file.exists(expr_file)){
  scores <- fread(expr_file)
  # 检查是否有其他线粒体基因的列
  available_c1 <- complex1_genes[complex1_genes %in% names(scores)]
  message("Complex I genes in score file: ", paste(available_c1, collapse = ", "))
  
  if(length(available_c1) > 1){
    # 计算每个基因与Ferroptosis Defense的相关
    ferro_def <- scores$Ferroptosis_Defense
    results <- data.frame()
    for(g in available_c1){
      if(g %in% names(scores)){
        r <- cor(scores[[g]], ferro_def, method = "spearman", use = "pairwise.complete.obs")
        results <- rbind(results, data.frame(
          Gene = g,
          Is_NDUFB7 = g == "NDUFB7",
          Spearman_Rho = round(r, 3),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    if(nrow(results) > 0){
      results <- results[order(abs(Spearman_Rho), decreasing = TRUE), ]
      fwrite(results, file.path(outdir, "V167C_complex1_ferroptosis_correlation.csv"))
      message("\n=== Complex I基因与Ferroptosis Defense相关 ===")
      print(results)
      
      ndufb7_r <- results$Spearman_Rho[results$Gene == "NDUFB7"]
      others_r <- results$Spearman_Rho[results$Gene != "NDUFB7"]
      message("\nNDUFB7 rho: ", ndufb7_r)
      message("Other Complex I median rho: ", round(median(abs(others_r), na.rm = TRUE), 3))
      
      if(abs(ndufb7_r) > max(abs(others_r), na.rm = TRUE)){
        message("[PASS] NDUFB7 has the STRONGEST correlation with ferroptosis defense among Complex I")
      } else {
        message("[WARN] Another Complex I gene has stronger correlation — NDUFB7 may not be specific")
      }
    }
  }
} else {
  message("[FAIL] Score file not found")
}

message("\n[DONE] V167C: ", outdir)
