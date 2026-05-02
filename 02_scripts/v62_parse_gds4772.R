#!/usr/bin/env Rscript
# V62: GDS4772 暴力解析 — 不依赖 GPL6244 注释包
# 直接从 GDS soft 文件提取 NDUFB7 探针 + 统计值

library(data.table)
library(GEOquery)

PROJECT <- Sys.getenv("HOME", "/home/y411869")
GDS_PATH <- file.path(PROJECT, "Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GDS4772/GDS4772.soft")
OUT_DIR <- file.path(PROJECT, "Projects/NDUFB7_HF_2026_04_20/03_results/02_tables")

cat("=== V62 GDS4772 暴力解析 ===\n")
cat("文件:", GDS_PATH, "\n")

# 方法1: 尝试用 GEOquery 直接解析 (如果可用)
if (requireNamespace("GEOquery", quietly = TRUE)) {
    cat("\n[方法1] 尝试 GEOquery::getGEO(filename=) ...\n")
    tryCatch({
        gds <- getGEO(filename = GDS_PATH, getGPL = FALSE)
        cat("✅ GEOquery 解析成功\n")
        
        # 提取表达矩阵
        expr <- Table(gds)
        cat("表达矩阵维度:", dim(expr), "\n")
        
        # 搜索 NDUFB7
        # 列名可能是 ID_REF, IDENTIFIER, Gene Symbol 等
        ndufb7_rows <- expr[grepl("NDUFB7", expr[,1], ignore.case = TRUE) | 
                           grepl("NDUFB7", expr[,2], ignore.case = TRUE), ]
        cat("NDUFB7 匹配行数:", nrow(ndufb7_rows), "\n")
        print(head(ndufb7_rows))
        
        # 保存
        fwrite(ndufb7_rows, file.path(OUT_DIR, "GDS4772_NDUFB7_raw_extract.csv"))
        cat("✅ 已保存: GDS4772_NDUFB7_raw_extract.csv\n")
    }, error = function(e) {
        cat("❌ GEOquery 解析失败:", conditionMessage(e), "\n")
    })
}

# 方法2: 暴力文本解析 (如果 GEOquery 失败)
cat("\n[方法2] 暴力文本解析...\n")
lines <- readLines(GDS_PATH, warn = FALSE)

# 找表达矩阵起始
start_idx <- grep("^!dataset_table_begin", lines)
end_idx <- grep("^!dataset_table_end", lines)

if (length(start_idx) > 0 && length(end_idx) > 0) {
    start_idx <- tail(start_idx, 1)
    end_idx <- tail(end_idx, 1)
    cat("矩阵范围: 行", start_idx+1, "到", end_idx-1, "\n")
    
    # 提取矩阵内容
    matrix_lines <- lines[(start_idx+1):(end_idx-1)]
    
    # 第一行是表头
    header <- strsplit(matrix_lines[1], "\t")[[1]]
    cat("表头列数:", length(header), "\n")
    cat("样本列:", length(header)-1, "\n")
    
    # 搜索 NDUFB7
    ndufb7_lines <- matrix_lines[grepl("NDUFB7", matrix_lines, ignore.case = TRUE)]
    cat("NDUFB7 匹配行数:", length(ndufb7_lines), "\n")
    
    if (length(ndufb7_lines) > 0) {
        # 解析为数据框
        vals <- strsplit(ndufb7_lines[1], "\t")[[1]]
        probe_id <- vals[1]
        expr_vals <- as.numeric(vals[-1])
        
        cat("\n探针ID:", probe_id, "\n")
        cat("表达值统计:\n")
        cat("  n =", length(expr_vals), "\n")
        cat("  mean =", round(mean(expr_vals, na.rm = TRUE), 4), "\n")
        cat("  sd =", round(sd(expr_vals, na.rm = TRUE), 4), "\n")
        cat("  median =", round(median(expr_vals, na.rm = TRUE), 4), "\n")
        cat("  range =", round(range(expr_vals, na.rm = TRUE), 4), "\n")
        
        # 样本分组信息 (从 soft 文件提取)
        # 找样本标题和描述
        title_lines <- grep("^\\^SAMPLE = ", lines, value = TRUE)
        desc_lines <- grep("!Sample_description", lines, value = TRUE)
        cat("\n样本数:", length(title_lines), "\n")
        
        # 简单分组: 搜索 DCM/NF 关键词
        dcm_idx <- grep("DCM|dilated|cardiomyopathy", lines[(start_idx+1):end_idx], ignore.case = TRUE)
        nf_idx <- grep("non-failing|normal|control|healthy", lines[(start_idx+1):end_idx], ignore.case = TRUE)
        
        # 保存结果
        result <- data.frame(
            probe_id = probe_id,
            mean_expr = mean(expr_vals, na.rm = TRUE),
            sd_expr = sd(expr_vals, na.rm = TRUE),
            median_expr = median(expr_vals, na.rm = TRUE),
            n_samples = length(expr_vals),
            note = "GDS4772_direct_extraction_no_GPL"
        )
        fwrite(result, file.path(OUT_DIR, "GDS4772_NDUFB7_stats.csv"))
        cat("\n✅ 统计值已保存: GDS4772_NDUFB7_stats.csv\n")
    } else {
        cat("❌ 未在表达矩阵中找到 NDUFB7\n")
        cat("建议: 检查探针ID是否为 ILMN_ 格式或其他命名\n")
    }
} else {
    cat("❌ 未找到 dataset_table_begin/end 标记\n")
}

cat("\n=== 解析完成 ===\n")
