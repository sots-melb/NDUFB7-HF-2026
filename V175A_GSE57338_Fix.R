suppressPackageStartupMessages(library(data.table))
library(dplyr)
PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V175A: GSE57338 修复")
message("========================================")

# 1. 读取表达矩阵
expr_file <- "01_data/01_raw_geo/GSE57338/GSE57338_series_matrix_expr.txt"
if(!file.exists(expr_file)) {
  message("[FAIL] expr file not found")
  quit(status=1)
}

# 先诊断格式
df_preview <- fread(expr_file, nrows=3, header=TRUE, sep="\t", quote="")
message("[DIAG] expr format: ", nrow(df_preview), " x ", ncol(df_preview))
message("  First col: ", names(df_preview)[1])
message("  First ID: ", df_preview[[1]][1])

# 完整读取
expr <- fread(expr_file, header=TRUE, sep="\t", quote="", fill=TRUE)
probe_ids <- expr[[1]]
expr_mat <- as.matrix(expr[, -1, with=FALSE])
rownames(expr_mat) <- probe_ids
message("[PASS] Expression matrix: ", nrow(expr_mat), " x ", ncol(expr_mat))

# 2. 解析GPL11532（用comment.char跳过#注释）
gpl_file <- "01_data/01_raw_geo/GPL11532/GPL11532.annot.gz"
if(!file.exists(gpl_file)) {
  message("[FAIL] GPL11532 not found")
  quit(status=1)
}

gpl <- fread(gpl_file, comment.char="#", header=TRUE, sep="\t", quote="", fill=TRUE)
message("[PASS] GPL11532 parsed: ", nrow(gpl), " x ", ncol(gpl))
message("  Columns: ", paste(names(gpl)[1:min(5, ncol(gpl))], collapse=", "))

# 找ID列和Gene Symbol列
id_col <- names(gpl)[1]
gene_candidates <- grep("Gene|Symbol|gene_symbol|GENE_SYMBOL", names(gpl), value=TRUE, ignore.case=TRUE)
gene_col <- ifelse(length(gene_candidates)>0, gene_candidates[1], names(gpl)[2])

message("  ID col: ", id_col)
message("  Gene col: ", gene_col)

mapping <- gpl[, c(id_col, gene_col), with=FALSE]
setnames(mapping, c("ProbeID", "GeneSymbol"))
mapping <- mapping[!is.na(GeneSymbol) & GeneSymbol != "" & GeneSymbol != "---" & GeneSymbol != "null"]
mapping <- mapping[!duplicated(ProbeID)]

message("[PASS] Mapping table: ", nrow(mapping), " entries")

# 3. 映射探针→基因
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

# 4. 检查NDUFB7
if("NDUFB7" %in% rownames(gene_mat)) {
  vals <- gene_mat["NDUFB7",]
  message("[FOUND] NDUFB7. Range: ", round(min(vals),3), " - ", round(max(vals),3), ", mean ", round(mean(vals),3))
  saveRDS(gene_mat, "03_results/V175_Fixed/GSE57338_gene_level_V175.rds")
  write.csv(data.frame(Gene=rownames(gene_mat), gene_mat, check.names=FALSE),
            "03_results/V175_Fixed/GSE57338_gene_level_V175.csv", row.names=FALSE)
  message("[DONE] Saved to 03_results/V175_Fixed/")
} else {
  message("[MISS] NDUFB7 not found")
  ndufb7_like <- grep("NDUFB|NDUF", rownames(gene_mat), value=TRUE)
  if(length(ndufb7_like) > 0) message("  Similar: ", paste(head(ndufb7_like, 5), collapse=", "))
}
