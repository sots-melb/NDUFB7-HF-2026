suppressMessages(library(GEOquery))
mat <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE154170/GSE154170_series_matrix.txt.gz"
if(file.exists(mat)) {
  gse <- getGEO(filename=mat, GSEMatrix=TRUE, getGPL=FALSE, AnnotGPL=FALSE)
  pdata <- pData(gse)
  message("▶ GSE154170 临床特征列:")
  surv_cols <- grep("characteristics_ch1", colnames(pdata), value=TRUE)
  if(length(surv_cols)>0) {
    for(c in surv_cols) {
      message("  列 ", c, ": "); print(head(pdata[[c]], 5))
    }
  }
  message("▶ 全部分组信息:")
  print(table(pdata$source_name_ch1, useNA="ifany"))
  message("▶ 标题信息:")
  print(head(pdata$title, 10))
} else {
  message("❌ series matrix不存在")
}
