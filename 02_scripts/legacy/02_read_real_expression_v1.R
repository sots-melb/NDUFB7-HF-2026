# ==============================================================================
# 脚本: 02_read_real_expression_v1.R
# 功能：读取GEO下载的真实表达数据（RAW.tar / CSV / Series Matrix）
# 针对：GSE141910_RAW.tar, GSE168742 CSV, 其他Series Matrix
# ==============================================================================

PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

library(tidyverse)

# ------------------- GSE141910：从RAW.tar读取 -------------------
# GSE141910是RNA-seq，RAW.tar内是CSV格式的counts/FPKM文件
cat("\n========== [处理] GSE141910 RAW数据 ==========\n")

raw_tar <- file.path(PROJECT_ROOT, "01_data/01_raw_geo/GSE141910/GSE141910_RAW.tar")
if (file.exists(raw_tar)) {
  # 列出tar内文件
  tar_files <- untar(raw_tar, list = TRUE)
  cat("[信息] RAW.tar内含文件数:", length(tar_files), "\n")
  cat("[信息] 前5个文件:\n")
  print(head(tar_files, 5))
  
  # 读取第一个文件查看格式（通常是CSV）
  # 注意：untar到临时目录读取，不在原目录解压
  temp_dir <- file.path(PROJECT_ROOT, "06_tmp/gse141910_temp")
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 只解压第一个文件查看结构
  untar(raw_tar, files = tar_files[1], exdir = temp_dir)
  sample1_file <- list.files(temp_dir, full.names = TRUE)[1]
  
  if (!is.na(sample1_file)) {
    # 自动检测分隔符
    first_lines <- readLines(sample1_file, n = 3)
    cat("[信息] 文件前3行:\n")
    print(first_lines)
    
    # 根据格式读取（可能是tab或逗号分隔）
    if (grepl("\t", first_lines[1])) {
      df1 <- read.delim(sample1_file, header = TRUE, row.names = 1)
    } else {
      df1 <- read.csv(sample1_file, header = TRUE, row.names = 1)
    }
    
    cat("[成功] 示例样本:", ncol(df1), "列（基因）x", nrow(df1), "行\n")
    cat("[成功] NDUFB7在此样本中的表达量:", 
        as.numeric(df1["NDUFB7", 1]), "\n")
  }
  
  # 清理临时文件（节省空间）
  unlink(temp_dir, recursive = TRUE)
  
} else {
  cat("[错误] GSE141910_RAW.tar 不存在，请先执行下载脚本\n")
}

# ------------------- GSE168742：单细胞CSV读取 -------------------
cat("\n========== [处理] GSE168742 单细胞数据 ==========\n")

sc_files <- c(
  HF = file.path(PROJECT_ROOT, "01_data/01_raw_geo/GSE168742/GSE168742_human_HF_CM.csv.gz"),
  Control = file.path(PROJECT_ROOT, "01_data/01_raw_geo/GSE168742/GSE168742_human_control_CM.csv.gz")
)

for (name in names(sc_files)) {
  f <- sc_files[[name]]
  if (file.exists(f)) {
    cat("[读取]", name, "...\n")
    # 使用data.table快速读取大CSV
    if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
    library(data.table)
    
    dt <- fread(f, nrows = 5)  # 先读5行看结构
    cat("[成功]", name, "结构:", ncol(dt), "列 x", nrow(dt), "行\n")
    cat("[成功] 列名前5个:", paste(head(names(dt), 5), collapse = ", "), "\n")
    
    # 检查NDUFB7
    if ("NDUFB7" %in% names(dt)) {
      cat("[QC-通过] NDUFB7在列名中存在！\n")
    } else {
      cat("[QC-警告] NDUFB7不在列名中，可能需要行名转换\n")
    }
  } else {
    cat("[跳过]", name, "文件不存在:", f, "\n")
  }
}

cat("\n[完成] 真实数据读取脚本结束\n")
cat("[下一步] 根据GSE141910_RAW.tar内文件格式，决定是否批量合并为表达矩阵\n")
