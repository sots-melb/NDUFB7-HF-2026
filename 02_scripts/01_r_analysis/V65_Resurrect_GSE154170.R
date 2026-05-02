suppressMessages(library(GEOquery))
message("▶ 开始解析 GSE154170 带有平台后缀的矩阵...")
base_dir <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE154170"

# 动态寻找所有矩阵文件 (GPL24676, GPL24247 等)
files <- list.files(path = base_dir, pattern = "series_matrix.txt.gz$", full.names = TRUE)

if(length(files) > 0) {
    # 优先解析第一个矩阵以获取临床表型
    gse <- getGEO(filename=files[1], GSEMatrix=TRUE)
    if(is.list(gse)) gse <- gse[[1]]
    
    pdata <- pData(gse)
    
    # 动态判定分组：根据 title 或 characteristics_ch1 识别 ICM, DCM 或 Control
    # 这里的 grep 条件需要尽可能宽泛以捕获不同写法
    pdata$Disease_Group <- "Unknown"
    pdata$Disease_Group[grepl("RSB|Control|Healthy|Non-failing", pdata$title, ignore.case=TRUE)] <- "Control"
    pdata$Disease_Group[grepl("ICM|Ischemic", pdata$title, ignore.case=TRUE)] <- "ICM"
    pdata$Disease_Group[grepl("DCM|Dilated", pdata$title, ignore.case=TRUE)] <- "DCM"
    
    message("\n✅ 分组判定完成 (GSE154170)。当前各组样本数量：")
    print(table(pdata$Disease_Group))
    
    out_csv <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE154170_Resurrected_Group.csv"
    write.csv(pdata, out_csv)
    message(paste("🎉 临床表型提取大获成功，已保存至:", out_csv))
} else {
    message("❌ 转移失败，未在目标文件夹找到矩阵。")
}
