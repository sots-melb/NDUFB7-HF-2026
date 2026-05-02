# ==============================================================================
# 脚本: 04_verify_ndufb7_real_data.R
# 功能: 读取GSE141910 RAW.tar + GSE168742 CSV，验证NDUFB7存在性
# 日期: 2026-04-20
# ==============================================================================

PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

# ------------------- GSE168742: 单细胞CSV（最优先，已确认真实） -------------------
cat("\n========== [读取] GSE168742 单细胞数据 ==========\n")
library(data.table)

# 读取HF样本（17M）
hf_file <- "01_data/01_raw_geo/GSE168742/GSE168742_human_HF_CM.csv.gz"
if (file.exists(hf_file)) {
  hf_dt <- fread(hf_file, nrows = 10)  # 先读10行看结构
  cat("HF文件结构: ", ncol(hf_dt), "列 x ", nrow(hf_dt), "行 (仅前10行)\n")
  cat("列名前10个: ", paste(head(names(hf_dt), 10), collapse = ", "), "\n")
  cat("行名前10个: ", paste(head(rownames(hf_dt), 10), collapse = ", "), "\n")
  
  # NDUFB7检测（可能在行名或列名）
  if ("NDUFB7" %in% names(hf_dt)) {
    cat("✅ NDUFB7在列名中找到！\n")
  } else if ("NDUFB7" %in% rownames(hf_dt)) {
    cat("✅ NDUFB7在行名中找到！\n")
  } else {
    cat("⚠️ NDUFB7未直接找到，尝试模糊匹配...\n")
    print(grep("NDUFB7", names(hf_dt), value = TRUE, ignore.case = TRUE)[1:5])
  }
  
  # 读取完整数据（如果结构确认正确）
  cat("\n[读取] 正在读取完整HF数据（17MB，可能需要10-30秒）...\n")
  hf_full <- fread(hf_file)
  cat("[完成] 完整HF: ", ncol(hf_full), "列 x ", nrow(hf_full), "行\n")
  saveRDS(hf_full, "01_data/01_raw_geo/GSE168742/GSE168742_human_HF_CM.rds")
  cat("[保存] GSE168742_human_HF_CM.rds\n")
} else {
  cat("❌ HF文件不存在\n")
}

# 读取Control样本（1.3M）
ctrl_file <- "01_data/01_raw_geo/GSE168742/GSE168742_human_control_CM.csv.gz"
if (file.exists(ctrl_file)) {
  ctrl_full <- fread(ctrl_file)
  cat("\n[完成] 完整Control: ", ncol(ctrl_full), "列 x ", nrow(ctrl_full), "行\n")
  saveRDS(ctrl_full, "01_data/01_raw_geo/GSE168742/GSE168742_human_control_CM.rds")
  cat("[保存] GSE168742_human_control_CM.rds\n")
}

# ------------------- GSE141910: RAW.tar（批量读取） -------------------
cat("\n========== [读取] GSE141910 RAW.tar ==========\n")
raw_tar <- "01_data/01_raw_geo/GSE141910/GSE141910_RAW.tar"
if (file.exists(raw_tar)) {
  # 列出tar内文件
  tar_list <- untar(raw_tar, list = TRUE)
  cat("RAW.tar内含文件数: ", length(tar_list), "\n")
  cat("前5个文件:\n")
  print(head(tar_list, 5))
  
  # 读取第一个文件确定格式
  temp_dir <- "06_tmp/gse141910_probe"
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  untar(raw_tar, files = tar_list[1], exdir = temp_dir)
  probe_file <- list.files(temp_dir, full.names = TRUE)[1]
  
  if (!is.na(probe_file)) {
    probe_lines <- readLines(probe_file, n = 3)
    cat("\n第一个文件前3行:\n")
    print(probe_lines)
    
    # 根据分隔符读取
    if (grepl("\t", probe_lines[1])) {
      df_probe <- read.delim(probe_file, header = TRUE, row.names = 1)
    } else {
      df_probe <- read.csv(probe_file, header = TRUE, row.names = 1)
    }
    cat("\n示例样本维度: ", ncol(df_probe), "列 x ", nrow(df_probe), "行\n")
