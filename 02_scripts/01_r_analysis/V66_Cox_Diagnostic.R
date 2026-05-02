suppressMessages(library(GEOquery))
matrix_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz"
gse <- getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE, AnnotGPL = FALSE)
pdata <- pData(gse)

message("==================================================")
message("▶ GSE59867 临床特征诊断面板")
message("==================================================")
# 寻找可能包含生存信息的列
surv_cols <- grep("characteristics_ch1", colnames(pdata), value=TRUE)
if(length(surv_cols) > 0) {
    message("找到以下特征列的前 3 行数据：")
    print(pdata[1:3, surv_cols, drop=FALSE])
} else {
    message("未找到 characteristics_ch1，打印所有列名：")
    print(colnames(pdata))
}
