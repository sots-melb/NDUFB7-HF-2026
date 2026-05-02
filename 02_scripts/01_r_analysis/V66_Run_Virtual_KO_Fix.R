suppressMessages(library(Seurat))
suppressMessages(library(scTenifoldKnk))

rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS"
obj <- readRDS(rds_path)

set.seed(42)
keep <- sample(Cells(obj), min(500, ncol(obj)))
obj_sub <- subset(obj, cells = keep)

# 核心修复：兼容 Seurat V4 和 V5 的 counts 提取方式
counts <- tryCatch({
    GetAssayData(obj_sub, slot = "counts") # V4
}, error = function(e) {
    GetAssayData(obj_sub, layer = "counts") # V5
})

# 转换为标准普通矩阵格式以防止 scTenifoldKnk 报错
counts <- as.matrix(counts)
counts <- counts[rowSums(counts) > 0, ]

message("▶ 矩阵提取成功！开始构建网络并进行 NDUFB7 敲除 (可能需要10-20分钟，请不要关闭终端)...")
ko_res <- scTenifoldKnk(countMatrix = counts, gKO = "NDUFB7", qc = FALSE)

dr_table <- ko_res$diffRegulation
dr_table <- dr_table[order(dr_table$p.adj), ]

out_file <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/NDUFB7_Virtual_KO_Predictions.csv"
write.csv(dr_table, out_file, row.names = FALSE)

message("🎉 虚拟 KO 分析大功告成！")
message("前 10 个受 NDUFB7 敲除显著影响的下游基因：")
print(head(dr_table, 10))
