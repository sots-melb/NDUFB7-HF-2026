suppressMessages(library(GEOquery))
for(gse_id in c("GSE57338", "GSE59867")) {
  message("\n▶ 检查 ", gse_id, " ...")
  mat <- list.files(paste0("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/", gse_id), pattern="series_matrix.txt.gz", full.names=TRUE)
  if(length(mat)>0) {
    gse <- getGEO(filename=mat[1], GSEMatrix=TRUE, getGPL=FALSE, AnnotGPL=FALSE)
    pdata <- pData(gse)
    cols <- colnames(pdata)
    surv_cols <- cols[grepl("time|follow|survival|day|month|year|status|event|death|outcome|endpoint|transplant", cols, ignore.case=TRUE)]
    if(length(surv_cols)>0) {
      message("  ✅ 发现潜在生存列: ", paste(surv_cols, collapse=", "))
      for(c in surv_cols) message("     - ", c, ": ", paste(head(pdata[[c]], 3), collapse=" | "))
    } else {
      message("  ❌ 未发现明显生存时间记录")
    }
  } else {
    message("  ❌ 本地未找到 series matrix")
  }
}
