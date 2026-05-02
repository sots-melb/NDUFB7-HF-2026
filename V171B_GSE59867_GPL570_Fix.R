suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V171B: GSE59867 GPL570探针修复")
message("========================================")

# 1. 解析series matrix获取表达值
sm_file <- "01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz"
if(!file.exists(sm_file)) {
  message("[FAIL] Series matrix not found")
  quit(status=1)
}

# 读取并定位表达矩阵
lines <- readLines(sm_file, n=100)
matrix_start <- grep("!series_matrix_table_begin", lines)
if(length(matrix_start) == 0) {
  # 尝试从压缩文件读取
  lines <- readLines(gzfile(sm_file), n=200)
  matrix_start <- grep("!series_matrix_table_begin", lines)
}

message("[DIAG] Matrix start line: ", matrix_start[1])

# 使用data.table fread（自动跳过注释行）
sm <- fread(sm_file, skip="!series_matrix_table_begin", header=TRUE, fill=TRUE)
message("[PASS] Loaded series matrix: ", nrow(sm), " x ", ncol(sm))

# 确定ID列和样本列
id_col <- names(sm)[1]
sample_cols <- names(sm)[2:ncol(sm)]

message("  ID column: ", id_col)
message("  Samples: ", length(sample_cols))

# 提取探针ID和表达矩阵
probe_ids <- sm[[id_col]]
expr_mat <- as.matrix(sm[, -1, with=FALSE])
rownames(expr_mat) <- probe_ids

message("[PASS] Expression matrix: ", nrow(expr_mat), " x ", ncol(expr_mat))
message("  Sample probes: ", paste(head(probe_ids, 3), collapse=", "), "...")

# 2. 查找GPL570注释
gpl_files <- c(
  "01_data/02_gpl_annotation/GPL570.annot.gz",
  "01_data/02_gpl_annotation/GPL570.annot",
  "Downloads/GPL570.annot.gz",
  "Downloads/GPL570.annot"
)

gpl_file <- NULL
for(f in gpl_files) {
  if(file.exists(f)) {
    gpl_file <- f
    message("[FOUND] GPL570: ", f)
    break
  }
}

if(is.null(gpl_file)) {
  message("[FAIL] GPL570.annot.gz not found. Need manual download.")
  message("  URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570")
  message("  Click: 'View full table' → 'Save all results to file' → download GPL570.annot.gz")
  
  # 保存探针级供后续
  saveRDS(expr_mat, "03_results/V171_Fixed/GSE59867_probe_level.rds")
  write.csv(data.frame(ProbeID=probe_ids), "03_results/V171_Fixed/GSE59867_probe_list.csv", row.names=FALSE)
  quit(status=1)
}

# 3. 解析GPL570
message("[STEP 3] Parsing GPL570...")
gpl <- fread(gpl_file, skip="#ID", header=TRUE, fill=TRUE, quote="")

id_col_gpl <- names(gpl)[1]
gene_col <- grep("Gene|Symbol|GENE_SYMBOL|gene_symbol", names(gpl), value=TRUE)[1]

if(is.na(gene_col)) {
  message("[DIAG] Available columns: ", paste(names(gpl)[1:min(10,ncol(gpl))], collapse=", "))
  gene_col <- names(gpl)[which(grepl("Gene|Symbol", names(gpl), ignore.case=TRUE))[1]]
}

message("  GPL ID col: ", id_col_gpl)
message("  GPL Gene col: ", gene_col)

mapping <- gpl[, c(id_col_gpl, gene_col), with=FALSE]
setnames(mapping, c("ProbeID", "GeneSymbol"))
mapping <- mapping[!is.na(GeneSymbol) & GeneSymbol != "" & GeneSymbol != "---"]
mapping <- mapping[!duplicated(ProbeID)]

message("[PASS] GPL570 mapping: ", nrow(mapping), " entries")

# 4. 映射到基因
matched <- probe_ids %in% mapping$ProbeID
message("[STEP 4] Matching: ", sum(matched), " / ", length(probe_ids))

expr_mapped <- expr_mat[matched, , drop=FALSE]
gene_symbols <- mapping$GeneSymbol[match(rownames(expr_mapped), mapping$ProbeID)]

gene_df <- data.frame(GeneSymbol=gene_symbols, expr_mapped, check.names=FALSE)
gene_unique <- gene_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), median), .groups="drop")

gene_mat <- as.matrix(gene_unique[,-1])
rownames(gene_mat) <- gene_unique$GeneSymbol

message("[PASS] Gene matrix: ", nrow(gene_mat), " x ", ncol(gene_mat))

# 5. 提取目标基因
target_genes <- c("NDUFB7", "ACSL4", "GPX4", "SLC7A11", "FTH1", "FTL", 
                  "NFE2L2", "KEAP1", "NOX1", "NOX2", "NOX4", "SOD1", "SOD2",
                  "CAT", "PRDX1", "TXNRD1", "GCLC", "GCLM", "HMOX1")

found <- target_genes %in% rownames(gene_mat)
message("[STEP 5] Target genes found: ", sum(found), "/", length(target_genes))
message("  Found: ", paste(target_genes[found], collapse=", "))
if(sum(!found) > 0) message("  Missing: ", paste(target_genes[!found], collapse=", "))

# 保存基因级矩阵
saveRDS(gene_mat, "03_results/V171_Fixed/GSE59867_gene_level_FIXED.rds")

# 6. 如果找到NDUFB7，做基础统计
if("NDUFB7" %in% rownames(gene_mat)) {
  ndufb7_vals <- gene_mat["NDUFB7", ]
  message("[FOUND] NDUFB7 range: ", round(min(ndufb7_vals), 3), " - ", round(max(ndufb7_vals), 3))
}

# 保存提取结果
extracted <- gene_mat[target_genes[found], , drop=FALSE]
if(nrow(extracted) > 0) {
  write.csv(data.frame(Gene=rownames(extracted), extracted, check.names=FALSE),
            "03_results/V171_Fixed/GSE59867_target_genes.csv", row.names=FALSE)
  message("[DONE] Saved extracted genes")
}
