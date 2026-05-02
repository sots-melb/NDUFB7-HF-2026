if(!require("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
library(GEOquery)

message("开始解析 GSE59867 Series Matrix 获取预后随访信息...")
tryCatch({
    # 获取 GSE59867 临床矩阵
    gse <- getGEO("GSE59867", GSEMatrix = TRUE, getGPL = FALSE)
    pdata <- pData(gse[[1]])
    
    # 提取需要的变量 (Time, Event, Age, Sex 等)
    # 注意: 实际列名需要根据 pdata 的 names 动态 grep
    write.csv(pdata, "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE59867_Clinical_Raw.csv")
    message("✅ GSE59867 临床表型提取成功，准备就绪用于 Cox 回归。")
}, error = function(e) {
    message("❌ GSE59867 提取失败，请检查网络或 FTP (E89教训): ", e$message)
})
