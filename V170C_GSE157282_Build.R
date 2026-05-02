suppressPackageStartupMessages(library(data.table))
PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V170C: GSE157282 DOX模型矩阵构建")
message("========================================")

tar_file <- "Downloads/GSE157282_RAW.tar"
if(!file.exists(tar_file)) tar_file <- "01_data/01_raw_geo/GSE157282/GSE157282_RAW.tar"

if(!file.exists(tar_file)) {
  message("[FAIL] GSE157282_RAW.tar not found")
  quit(save="no", status=1)
}

message("\n[1/3] Extracting tar...")
tmp_dir <- tempfile(pattern = "gse157282_")
dir.create(tmp_dir, showWarnings = FALSE)
system(paste("tar -xf", tar_file, "-C", tmp_dir, "--wildcards", "*.txt.gz"))

files <- list.files(tmp_dir, pattern = "GSM.*\\.txt\\.gz$", full.names = TRUE)
message("  Found ", length(files), " sample files")

if(length(files) == 0) {
  message("[FAIL] No txt.gz files in tar")
  quit(save="no", status=1)
}

# 读取并合并
message("\n[2/3] Merging expression matrix...")
first <- fread(files[1], header = FALSE, col.names = c("ENSG", "Symbol", "Chr", "Value"))
mat <- data.table(Symbol = first$Symbol, S1 = first$Value)

for(i in 2:length(files)) {
  tmp <- fread(files[i], header = FALSE, col.names = c("ENSG", "Symbol", "Chr", "Value"))
  mat <- merge(mat, tmp[, .(Symbol, Value)], by = "Symbol", all = TRUE)
  setnames(mat, "Value", paste0("S", i))
  if(i %% 3 == 0) message("  Processed ", i, "/", length(files))
}

message("\n  Final: ", nrow(mat), " genes x ", ncol(mat)-1, " samples")

# 保存
saveRDS(mat, "03_results/V170_Rescue/GSE157282_rescued_matrix.rds")
fwrite(mat, "03_results/V170_Rescue/GSE157282_rescued_matrix.csv")

# 提取NDUFB7
nduf <- mat[tolower(Symbol) == "ndufb7"]
if(nrow(nduf) > 0) {
  vals <- as.numeric(nduf[1, -1])
  message("\n[PASS] NDUFB7 found")
  message("  Values: ", paste(round(vals, 2), collapse = ", "))
  
  # 简单分组比较（前3=Control, 后3=DOX，假设文件按顺序）
  # 实际分组需要从series matrix解析，这里仅做初步
  ctrl_mean <- mean(vals[1:3], na.rm = TRUE)
  dox_mean <- mean(vals[10:12], na.rm = TRUE)  # 根据V168A，DOX是最后3个
  message("  Control mean: ", round(ctrl_mean, 3))
  message("  DOX mean: ", round(dox_mean, 3))
  message("  Fold change: ", round(dox_mean / ctrl_mean, 3))
} else {
  message("\n[WARN] NDUFB7 not found by symbol")
}

unlink(tmp_dir, recursive = TRUE)
message("\n[DONE] V170C")
