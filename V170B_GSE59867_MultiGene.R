suppressPackageStartupMessages(library(data.table))
PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V170B: GSE59867 多基因提取 (ACSL4/GPX4/ROS/Death)")
message("========================================")

# === 1. 解析series matrix ===
sm_file <- "01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz"
if(!file.exists(sm_file)) sm_file <- "Downloads/GSE59867_series_matrix.txt.gz"

if(!file.exists(sm_file)) {
  message("[FAIL] GSE59867 series matrix not found")
  quit(save="no", status=1)
}

message("\n[1/5] Parsing: ", sm_file)
all_lines <- readLines(sm_file, n = 20000)
start_line <- which(grepl("^!series_matrix_table_begin", all_lines))
if(length(start_line) == 0) start_line <- which(grepl("^ID_REF", all_lines))
end_line <- which(grepl("^!series_matrix_table_end", all_lines))
if(length(end_line) == 0) end_line <- length(all_lines) + 1

message("  Matrix start: line ", start_line)
message("  Matrix end: line ", end_line)

if(length(start_line) == 0) {
  message("[FAIL] Cannot locate matrix in series matrix")
  quit(save="no", status=1)
}

nrows <- end_line - start_line - 1
mat <- fread(sm_file, skip = start_line, nrows = nrows, header = TRUE)
message("[2/5] Loaded: ", nrow(mat), " probes x ", ncol(mat)-1, " samples")

id_col <- names(mat)[1]
probe_ids <- as.character(mat[[id_col]])
message("  ID column: ", id_col)
message("  First probes: ", paste(head(probe_ids, 5), collapse = ", "))

# === 2. 获取GPL570注释 ===
message("\n[3/5] Resolving GPL570 annotation...")

# 尝试方法1: hgu133plus2.db
annot <- NULL
tryCatch({
  suppressPackageStartupMessages(library(hgu133plus2.db))
  annot <- select(hgu133plus2.db, keys = probe_ids, 
                  columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
  message("  [PASS] hgu133plus2.db annotation: ", nrow(annot), " mappings")
}, error = function(e) {
  message("  [WARN] hgu133plus2.db unavailable: ", conditionMessage(e))
})

# 尝试方法2: GPL570.annot.gz（如果已下载）
if(is.null(annot)) {
  gpl_file <- "01_data/01_raw_geo/GPL570/GPL570.annot.gz"
  if(!file.exists(gpl_file)) gpl_file <- "Downloads/GPL570.annot.gz"
  if(file.exists(gpl_file)) {
    message("  [TRY] Parsing GPL570.annot.gz...")
    gpl <- fread(gpl_file, skip = "#ID", header = TRUE, quote = "")
    # 找探针ID和Symbol列
    probe_col <- grep("ID|Probe", names(gpl), value = TRUE, ignore.case = TRUE)[1]
    sym_col <- grep("Symbol|Gene|SYMBOL", names(gpl), value = TRUE, ignore.case = TRUE)[1]
    if(!is.na(probe_col) && !is.na(sym_col)) {
      annot <- gpl[, c(probe_col, sym_col), with = FALSE]
      names(annot) <- c("PROBEID", "SYMBOL")
      message("  [PASS] GPL570.annot.gz: ", nrow(annot), " mappings")
    }
  }
}

# 尝试方法3: 从series matrix的探针ID推断（如果ID已经是基因名）
if(is.null(annot) && !any(grepl("_at$", head(probe_ids, 100)))) {
  message("  [DIAG] Probe IDs do not look like Affy IDs — assuming gene symbols")
  mat$SYMBOL <- probe_ids
} else if(!is.null(annot)) {
  # 合并注释
  annot <- annot[!is.na(SYMBOL) & SYMBOL != "", ]
  annot <- annot[!duplicated(PROBEID)]
  mat$PROBEID <- probe_ids
  mat <- merge(mat, annot[, .(PROBEID, SYMBOL)], by = "PROBEID", all.x = TRUE)
  
  # 对重复基因取中位数
  mat_gene <- mat[!is.na(SYMBOL), ]
  if(nrow(mat_gene) > 0) {
    sample_cols <- setdiff(names(mat_gene), c("PROBEID", "SYMBOL", "GENENAME"))
    # 按SYMBOL聚合（中位数）
    mat_gene <- mat_gene[, lapply(.SD, median, na.rm = TRUE), by = SYMBOL, .SDcols = sample_cols]
    message("  [PASS] Aggregated to gene level: ", nrow(mat_gene), " genes")
    mat <- mat_gene
    id_col <- "SYMBOL"
  }
}

# === 3. 提取目标基因 ===
message("\n[4/5] Extracting target genes...")

target_genes <- c("NDUFB7", "ACSL4", "GPX4", "SLC7A11", "FTH1", "FTL", 
                  "LPCAT3", "LOX", "PTGS2", "NOX4", "SOD2", "CAT", "GPX1",
                  "BCL2", "BAX", "CASP3", "RIPK1", "RIPK3", "MLKL")

gene_ids <- as.character(mat[[id_col]])
found <- target_genes[tolower(target_genes) %in% tolower(gene_ids)]
message("  Found: ", length(found), "/", length(target_genes))
message("  ", paste(found, collapse = ", "))

if(length(found) > 0) {
  idx <- match(tolower(found), tolower(gene_ids))
  extracted <- mat[idx, ]
  fwrite(extracted, "03_results/V170_Rescue/GSE59867_multigene_extracted.csv")
  message("\n[PASS] Saved: 03_results/V170_Rescue/GSE59867_multigene_extracted.csv")
  
  # 计算ACSL4/GPX4比值（如果两者都存在）
  if(all(c("ACSL4", "GPX4") %in% found)) {
    acsl4_row <- which(tolower(gene_ids) == "acsl4")
    gpx4_row <- which(tolower(gene_ids) == "gpx4")
    acsl4_vals <- as.numeric(mat[acsl4_row, -1])
    gpx4_vals <- as.numeric(mat[gpx4_row, -1])
    ratio <- acsl4_vals / (gpx4_vals + 0.001)
    message("\n  ACSL4/GPX4 ratio range: ", round(min(ratio, na.rm = TRUE), 3), 
            " - ", round(max(ratio, na.rm = TRUE), 3))
  }
} else {
  message("\n[WARN] No target genes found")
}

# 保存完整基因级矩阵（供后续使用）
if(id_col == "SYMBOL" && nrow(mat) > 1000) {
  fwrite(mat, "03_results/V170_Rescue/GSE59867_gene_level_matrix.csv")
  message("[5/5] Full gene-level matrix saved")
}

message("\n[DONE] V170B")
