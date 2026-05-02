suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V175B: GSE59867 修复")
message("========================================")

# 1. 行号定位解析series matrix
sm_file <- "01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz"
lines <- readLines(gzfile(sm_file), n=300)
start_idx <- grep("!series_matrix_table_begin", lines, fixed=TRUE)[1]
end_idx <- grep("!series_matrix_table_end", lines, fixed=TRUE)[1]

if(is.na(start_idx)) {
  message("[FAIL] Cannot find matrix start")
  quit(status=1)
}

message("[DIAG] Matrix start line: ", start_idx, ", end: ", end_idx)

# 从start_idx+1开始读（表头行），用comment.char="!"跳过!series_matrix_table_end
sm <- fread(sm_file, skip=start_idx, header=TRUE, sep="\t", quote="", comment.char="!", fill=TRUE)
message("[PASS] Series matrix: ", nrow(sm), " x ", ncol(sm))

probe_ids <- sm[[1]]
expr_mat <- as.matrix(sm[, -1, with=FALSE])
rownames(expr_mat) <- probe_ids
message("[PASS] Expression: ", nrow(expr_mat), " x ", ncol(expr_mat))
message("  Sample probes: ", paste(head(probe_ids, 3), collapse=", "))

# 2. 解析GPL570（comment.char跳过#注释）
gpl_file <- "01_data/02_gpl_annotation/GPL570_full.txt"
if(!file.exists(gpl_file)) {
  message("[FAIL] GPL570 not found")
  quit(status=1)
}

gpl <- fread(gpl_file, comment.char="#", header=TRUE, sep="\t", quote="", fill=TRUE)
message("[PASS] GPL570 parsed: ", nrow(gpl), " x ", ncol(gpl))
message("  Columns: ", paste(names(gpl)[1:min(5, ncol(gpl))], collapse=", "))

id_col <- names(gpl)[1]
gene_candidates <- grep("Gene|Symbol|gene_symbol|GENE_SYMBOL", names(gpl), value=TRUE, ignore.case=TRUE)
gene_col <- ifelse(length(gene_candidates)>0, gene_candidates[1], names(gpl)[2])

message("  ID col: ", id_col)
message("  Gene col: ", gene_col)

mapping <- gpl[, c(id_col, gene_col), with=FALSE]
setnames(mapping, c("ProbeID", "GeneSymbol"))
mapping <- mapping[!is.na(GeneSymbol) & GeneSymbol != "" & GeneSymbol != "---" & GeneSymbol != "null"]
mapping <- mapping[!duplicated(ProbeID)]

message("[PASS] Mapping: ", nrow(mapping), " entries")

# 3. 映射
matched <- probe_ids %in% mapping$ProbeID
message("[STEP 3] Matched: ", sum(matched), " / ", length(probe_ids))

if(sum(matched) == 0) {
  message("[FAIL] No probes matched")
  message("  Expr IDs: ", paste(head(probe_ids, 5), collapse=", "))
  message("  Map IDs: ", paste(head(mapping$ProbeID, 5), collapse=", "))
  quit(status=1)
}

expr_mapped <- expr_mat[matched, , drop=FALSE]
gene_symbols <- mapping$GeneSymbol[match(rownames(expr_mapped), mapping$ProbeID)]

gene_df <- data.frame(GeneSymbol=gene_symbols, expr_mapped, check.names=FALSE)
gene_unique <- gene_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), median), .groups="drop")

gene_mat <- as.matrix(gene_unique[,-1])
rownames(gene_mat) <- gene_unique$GeneSymbol

message("[PASS] Gene-level matrix: ", nrow(gene_mat), " x ", ncol(gene_mat))

# 4. 提取目标基因
targets <- c("NDUFB7", "ACSL4", "GPX4", "SLC7A11", "FTH1", "FTL", 
             "NFE2L2", "KEAP1", "NOX1", "NOX2", "NOX4", "SOD1", "SOD2",
             "CAT", "PRDX1", "TXNRD1", "GCLC", "GCLM", "HMOX1",
             "TP53", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "PARP1",
             "MLKL", "RIPK3", "RIPK1", "PGAM5", "GSDMD", "GSDME", "IL1B")

found <- targets %in% rownames(gene_mat)
message("[STEP 4] Target genes found: ", sum(found), "/", length(targets))
message("  Found: ", paste(targets[found], collapse=", "))
if(sum(!found)>0) message("  Missing: ", paste(targets[!found], collapse=", "))

saveRDS(gene_mat, "03_results/V175_Fixed/GSE59867_gene_level_V175.rds")
if(sum(found)>0) {
  extracted <- gene_mat[targets[found], , drop=FALSE]
  write.csv(data.frame(Gene=rownames(extracted), extracted, check.names=FALSE),
            "03_results/V175_Fixed/GSE59867_target_genes_V175.csv", row.names=FALSE)
}

if("NDUFB7" %in% rownames(gene_mat)) {
  vals <- gene_mat["NDUFB7",]
  message("[FOUND] NDUFB7: range ", round(min(vals),3), " - ", round(max(vals),3), ", mean ", round(mean(vals),3))
}

message("[DONE] Saved to 03_results/V175_Fixed/")
