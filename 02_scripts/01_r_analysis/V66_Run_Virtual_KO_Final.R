suppressMessages(library(Seurat))
suppressMessages(library(scTenifoldKnk))

rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS"
obj <- readRDS(rds_path)

set.seed(42)
keep <- sample(Cells(obj), min(500, ncol(obj)))
obj_sub <- subset(obj, cells = keep)

counts <- tryCatch({
    GetAssayData(obj_sub, slot = "counts")
}, error = function(e) {
    GetAssayData(obj_sub, layer = "counts")
})

counts <- as.matrix(counts)
counts <- counts[rowSums(counts) > 0, ]

message("▶ 矩阵提取成功！开始构建网络并进行 NDUFB7 敲除 (可能需要10-20分钟)...")
# 核心修复：去除了引起歧义的 qc=FALSE 参数
ko_res <- scTenifoldKnk(countMatrix = counts, gKO = "NDUFB7")

dr_table <- ko_res$diffRegulation
dr_table <- dr_table[order(dr_table$p.adj), ]

out_file <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/NDUFB7_Virtual_KO_Predictions.csv"
write.csv(dr_table, out_file, row.names = FALSE)

message("🎉 虚拟 KO 分析大功告成！")
message("前 10 个受影响的下游基因：")
print(head(dr_table, 10))
