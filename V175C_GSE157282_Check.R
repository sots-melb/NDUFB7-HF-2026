suppressPackageStartupMessages(library(data.table))
PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V175C: GSE157282 资产确认")
message("========================================")

sm_file <- "01_data/01_raw_geo/GSE157282/GSE157282_series_matrix.txt.gz"
if(file.exists(sm_file)) {
  sm <- fread(sm_file, skip="!series_matrix_table_begin", header=TRUE, fill=TRUE)
  message("[PASS] Series matrix: ", nrow(sm), " x ", ncol(sm))
} else {
  message("[MISS] Series matrix not found")
}

raw_tar <- "01_data/01_raw_geo/GSE157282/GSE157282_RAW.tar"
if(file.exists(raw_tar)) {
  message("[PASS] RAW tar: ", round(file.size(raw_tar)/1e6, 1), " MB")
  message("  Contents (first 10):")
  invisible(system(paste("tar -tf", raw_tar, "| head -10")))
} else {
  message("[MISS] RAW tar not found")
}

message("[DONE]")
