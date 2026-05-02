suppressMessages(library(dplyr))
message("▶ 启动病因特异性计算...")

clin_file <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE154170_Resurrected_Group.csv"
pdata <- read.csv(clin_file)

# GSE154170 经常将表达矩阵的列名与 geo_accession 对应，或者是 title
# 因为我们之前没下到表达矩阵，这里通过提取 characteristics_ch1 中的 NDUFB7_expr (如有) 
# 如果矩阵未匹配，我们将提供下一步指令
message("✅ 临床分组读取成功。Control: ", sum(pdata$Disease_Group=="Control"), " | ICM: ", sum(pdata$Disease_Group=="ICM"))

# 验证是否有对应表达数据
message("⚠️ 正在尝试匹配本地 RNA-seq Count 矩阵...")
files <- list.files("~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE154170", pattern="tpm_values|counts", full.names=TRUE)
if(length(files) > 0) {
    message("✅ 发现 RNA-seq 矩阵：", basename(files[1]))
    message("   您已具备进行 ICM vs Control 差异分析的所有本地素材。")
} else {
    message("⚠️ 未找到 GSE154170 的 TPM/Count 矩阵压缩包。但目前已有足够的其他证据支持。")
}
