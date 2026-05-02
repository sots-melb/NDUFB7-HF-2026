library(GEOquery)
message("检查本地 GSE59867_series_matrix.txt.gz...")
local_file <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz"
if(file.exists(local_file) && file.info(local_file)$size > 1000000) {
    gse <- getGEO(filename=local_file, GSEMatrix=TRUE)
    pdata <- pData(gse)
    write.csv(pdata, "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE59867_Clinical_Raw.csv")
    message("✅ GSE59867 本地提取完成！")
} else {
    message("⚠️ GSE59867 仍在后台下载中，稍后 R 脚本会自动重试...")
}
