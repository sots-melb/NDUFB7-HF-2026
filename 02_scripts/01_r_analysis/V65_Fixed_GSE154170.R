message("▶ 开始动态检索 GSE154170 ...")
base_dir <- "~/Projects/NDUFB7_HF_2026_04_20"

# 1. 动态搜索文件，绕过硬编码路径
files <- list.files(path = base_dir, pattern = "GSE154170.*series_matrix.txt.gz$", recursive = TRUE, full.names = TRUE)

if(length(files) > 0) {
    message(paste("✅ 成功定位矩阵文件:", files[1]))
    
    suppressMessages(library(GEOquery))
    # 2. 载入矩阵
    gse <- getGEO(filename=files[1], GSEMatrix=TRUE)
    
    # 3. 修复 S4 class 报错 (如果 getGEO 返回的是 list，强制提取第一个元素)
    if(is.list(gse)) {
        gse <- gse[[1]]
    }
    
    pdata <- pData(gse)
    
    # 4. 基于 V64 策略：使用样本名/标题前缀进行暴力分组 (RSB = Control, 其他 = ICM)
    # 假设 title 列包含类似 "RSB_xxx" 或 "ICM_xxx" 的信息
    pdata$Disease_Group <- ifelse(grepl("RSB", pdata$title, ignore.case=TRUE), "Control", "ICM")
    
    message("✅ 分组判定完成。当前各组样本数量：")
    print(table(pdata$Disease_Group))
    
    out_csv <- file.path(base_dir, "03_results/02_tables/GSE154170_Fixed_Group.csv")
    write.csv(pdata, out_csv)
    message(paste("✅ 修复后的临床表型已保存至:", out_csv))
} else {
    message("❌ 全局检索未发现 GSE154170_series_matrix.txt.gz，请确认该数据集是否曾成功下载。")
}
