suppressMessages(library(GEOquery))
options(timeout = 600)

datasets <- c("GSE120895", "GSE141910")
base_dir <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo"

for(gse in datasets) {
  message("▶ 正在探查 ", gse, " 的补充文件...")
  supp_files <- tryCatch(getGEOSuppFiles(gse, makeDirectory=FALSE, baseDir=base_dir, fetch_files=FALSE), error=function(e) NULL)
  
  if(!is.null(supp_files) && nrow(supp_files) > 0) {
    message("  发现补充文件:")
    print(supp_files$fname)
    
    # 下载非 RAW 的汇总表格 (通常是 xlsx, csv, txt)
    target_files <- supp_files$fname[grepl("xlsx|csv|txt|xls", supp_files$fname, ignore.case=TRUE) & !grepl("RAW", supp_files$fname)]
    
    if(length(target_files) > 0) {
      message("  正在下载疑似临床数据表: ", paste(target_files, collapse=", "))
      tryCatch(getGEOSuppFiles(gse, makeDirectory=FALSE, baseDir=base_dir, filter_regex="xlsx|csv|txt|xls"), error=function(e) message("  下载出错"))
    } else {
      message("  未发现表格型补充文件。")
    }
  } else {
    message("  无补充文件。")
  }
}
