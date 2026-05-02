suppressPackageStartupMessages(library(data.table))
PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V174A: GSE57338 直接修复")
message("========================================")

annotated_csv <- "03_results/02_tables/GSE57338_expression_matrix_annotated.csv"
if(file.exists(annotated_csv)) {
  message("[TRY] Loading annotated.csv...")
  df <- fread(annotated_csv, header=TRUE)
  message("  Dimensions: ", nrow(df), " x ", ncol(df))
  message("  First column: ", names(df)[1])
  
  gene_ids <- df[[1]]
  expr_mat <- as.matrix(df[, -1, with=FALSE])
  rownames(expr_mat) <- gene_ids
  
  message("  [PASS] Gene matrix: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")
  
  if("NDUFB7" %in% rownames(expr_mat)) {
    ndufb7_vals <- expr_mat["NDUFB7", ]
    message("[FOUND] NDUFB7. Range: ", round(min(ndufb7_vals),3), " - ", round(max(ndufb7_vals),3))
    saveRDS(expr_mat, "03_results/V174_Fixed/GSE57338_gene_level_DIRECT.rds")
    message("[DONE] Saved to 03_results/V174_Fixed/")
  } else {
    message("[MISS] NDUFB7 not found in annotated.csv")
    message("  Similar: ", paste(grep("NDUF", rownames(expr_mat), value=TRUE)[1:5], collapse=", "))
  }
} else {
  message("[FAIL] annotated.csv not found")
}
