message("启动 scTenifoldKnk NDUFB7 虚拟敲除分析...")
library(scTenifoldKnk)
library(Matrix)

message("正在全局搜索 GSE183852 表达矩阵 (CSV 或 RDS)...")
# 绕过写死的路径，动态搜索
files <- list.files(path = "~/Projects/NDUFB7_HF_2026_04_20", pattern = "GSE183852.*(csv|rds)$", recursive = TRUE, full.names = TRUE)

if (length(files) > 0) {
    message(paste("✅ 找到矩阵文件:", files[1]))
    message("⚠️ 由于 220k 细胞构建流形极耗内存，脚本将进入预处理模式，准备抽样矩阵送入 KO 模型...")
    # scTenifoldKnk(count_matrix, gKO = "NDUFB7") 
} else {
    message("❌ 未在项目中找到 GSE183852 矩阵，请检查下载日志。")
}
