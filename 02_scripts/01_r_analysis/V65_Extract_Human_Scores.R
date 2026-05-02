suppressMessages(library(Seurat))
message("▶ 仅载入 RDS 获取 metadata ...")
rds_file <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_figures_tables/GSE183852_scored.RDS"

if(file.exists(rds_file)) {
    # 尝试读取对象
    obj <- readRDS(rds_file)
    meta <- obj@meta.data
    
    # 尝试提取 NDUFB7 表达量
    if("NDUFB7" %in% rownames(obj)) {
        meta$NDUFB7_expr <- as.numeric(GetAssayData(obj, layer = "data")["NDUFB7", ])
        out_csv <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_extracted_meta.csv"
        write.csv(meta, out_csv)
        message(paste("🎉 提取大获成功！已保存至:", out_csv))
        
        # 简单查看已有评分列
        message("目前已有的打分列:")
        print(grep("Score", colnames(meta), value=TRUE))
    } else {
        message("❌ 对象中无 NDUFB7。")
    }
} else {
    message("❌ 文件不存在。")
}
