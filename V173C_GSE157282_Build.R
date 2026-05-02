suppressPackageStartupMessages(library(data.table))
PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V173C: GSE157282 矩阵构建")
message("========================================")

sm_file <- "01_data/01_raw_geo/GSE157282/GSE157282_series_matrix.txt.gz"
raw_tar <- "01_data/01_raw_geo/GSE157282/GSE157282_RAW.tar"

if(file.exists(sm_file)) {
  message("[FOUND] Series matrix")
  sm <- fread(sm_file, skip="!series_matrix_table_begin", header=TRUE, fill=TRUE)
  message("  Dimensions: ", nrow(sm), " x ", ncol(sm))
}

if(file.exists(raw_tar)) {
  message("[FOUND] RAW tar: ", round(file.size(raw_tar)/1e6,1), " MB")
  message("  Listing first 20 files:")
  invisible(system(paste("tar -tf", raw_tar, "| head -20")))
}

message("[DONE] GSE157282 assets ready for matrix construction")
