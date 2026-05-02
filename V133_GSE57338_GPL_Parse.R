#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V133_GSE57338")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V133_FIX: GPL11532.annot.gz暴力解析")
message("========================================")

gpl_file <- "~/Downloads/GPL11532.annot.gz"
if (!file.exists(gpl_file)) {
  # 尝试其他路径
  gpl_file <- list.files("~/Downloads", pattern="GPL11532.*\\.annot\\.gz$", full.names=TRUE)[1]
}
if (is.na(gpl_file) || !file.exists(gpl_file)) {
  message("[FAIL] GPL11532.annot.gz not found"); quit(status=1)
}

message("[LOAD] ", basename(gpl_file))

# 方法1: gunzip + readLines解析
tmp <- tempfile()
system(paste("gunzip -c", shQuote(gpl_file), ">", tmp))
lines <- readLines(tmp, n = 1000)
unlink(tmp)

# 探测表头
header_line <- grep("^ID", lines)[1]
if (is.na(header_line)) header_line <- grep("^\\^ID", lines)[1]
if (is.na(header_line)) {
  # 尝试直接data.table读取
  message("[INFO] 尝试data.table直接读取...")
  dt <- tryCatch(fread(paste0("zcat < ", shQuote(gpl_file)), skip="ID"), error=function(e) NULL)
  if (is.null(dt)) {
    message("[FAIL] 无法解析.annot.gz格式")
    quit(status=1)
  }
} else {
  message("[PASS] Header found at line ", header_line)
  # 提取列名
  header <- strsplit(lines[header_line], "\t")[[1]]
  message("Columns: ", paste(head(header, 5), collapse=", "), "...")
  
  # 读取完整文件
  dt <- fread(paste0("zcat < ", shQuote(gpl_file)), skip = header_line - 1, header = TRUE)
}

message("[PASS] Annot table: ", nrow(dt), " rows × ", ncol(dt), " cols")
message("Column names: ", paste(head(colnames(dt), 10), collapse=", "))

# 寻找基因名列
gene_col_candidates <- grep("Gene Symbol|gene_symbol|Symbol|GENE_SYMBOL", colnames(dt), value=TRUE, ignore.case=TRUE)
if (length(gene_col_candidates) == 0) {
  # 尝试其他常见命名
  gene_col_candidates <- grep("symbol|gene|title", colnames(dt), value=TRUE, ignore.case=TRUE)[1]
}
if (length(gene_col_candidates) > 0) {
  gene_col <- gene_col_candidates[1]
  message("[FOUND] Gene column: ", gene_col)
  
  # 提取NDUFB7探针
  ndufb7_rows <- dt[grep("NDUFB7", dt[[gene_col]], ignore.case=TRUE)]
  message("[PASS] NDUFB7 probes found: ", nrow(ndufb7_rows))
  if (nrow(ndufb7_rows) > 0) {
    write.csv(ndufb7_rows, file.path(outdir, "V133_NDUFB7_probes.csv"), row.names=FALSE)
    message("Probe IDs: ", paste(head(ndufb7_rows[[1]], 5), collapse=", "))
  }
  
  # 保存完整映射表
  probe_col <- colnames(dt)[1]  # 通常第一列是探针ID
  map_table <- dt[, c(probe_col, gene_col), with=FALSE]
  colnames(map_table) <- c("probe_id", "gene_symbol")
  write.csv(map_table, file.path(outdir, "V133_probe_to_gene_mapping.csv"), row.names=FALSE)
  message("[DONE] Full mapping saved: ", nrow(map_table), " probes")
  
  # 统计
  mapped_count <- sum(!is.na(map_table$gene_symbol) & map_table$gene_symbol != "")
  message("Mapped genes: ", mapped_count, " / ", nrow(map_table), 
          " (", round(mapped_count/nrow(map_table)*100,1), "%)")
  
} else {
  message("[FAIL] Cannot identify gene symbol column")
  message("[INFO] Available columns: ", paste(colnames(dt), collapse=", "))
}

message("[DONE] V133_FIX: ", outdir)
