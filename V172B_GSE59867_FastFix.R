suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V172B: GSE59867 GPL570 快速修复")
message("========================================")

# 1. 从series matrix直接解析（已知可用）
sm_file <- "01_data/01_raw_geo/GSE59867/GSE59867_series_matrix.txt.gz"
if(!file.exists(sm_file)) {
  message("[FAIL] Series matrix not found")
  quit(status=1)
}

# 读取表达矩阵部分
sm <- fread(sm_file, skip="!series_matrix_table_begin", header=TRUE, fill=TRUE, quote="")
message("[PASS] Loaded: ", nrow(sm), " x ", ncol(sm))

probe_ids <- sm[[1]]
expr_mat <- as.matrix(sm[, -1, with=FALSE])
rownames(expr_mat) <- probe_ids
message("[PASS] Expression matrix: ", nrow(expr_mat), " x ", ncol(expr_mat))

# 2. 加载GPL570
gpl_file <- "01_data/02_gpl_annotation/GPL570.annot.gz"
if(!file.exists(gpl_file)) {
  message("[FAIL] GPL570.annot.gz not found")
  quit(status=1)
}

message("[STEP 2] Parsing GPL570...")
gpl <- fread(gpl_file, skip="#ID", header=TRUE, fill=TRUE, quote="")

# 自动识别列
id_col <- names(gpl)[1]
gene_candidates <- grep("Gene|Symbol|GENE_SYMBOL|gene_symbol|Gene Symbol", names(gpl), value=TRUE)
gene_col <- ifelse(length(gene_candidates)>0, gene_candidates[1], names(gpl)[2])

message("  ID col: ", id_col)
message("  Gene col: ", gene_col)

mapping <- gpl[, c(id_col, gene_col), with=FALSE]
setnames(mapping, c("ProbeID", "GeneSymbol"))
mapping <- mapping[!is.na(GeneSymbol) & GeneSymbol != "" & GeneSymbol != "---"]
mapping <- mapping[!duplicated(ProbeID)]

message("[PASS] Mapping table: ", nrow(mapping), " entries")

# 3. 映射
matched <- probe_ids %in% mapping$ProbeID
message("[STEP 3] Matched: ", sum(matched), " / ", length(probe_ids))

expr_mapped <- expr_mat[matched, , drop=FALSE]
gene_symbols <- mapping$GeneSymbol[match(rownames(expr_mapped), mapping$ProbeID)]

gene_df <- data.frame(GeneSymbol=gene_symbols, expr_mapped, check.names=FALSE)
gene_unique <- gene_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), median), .groups="drop")

gene_mat <- as.matrix(gene_unique[,-1])
rownames(gene_mat) <- gene_unique$GeneSymbol

message("[PASS] Gene matrix: ", nrow(gene_mat), " x ", ncol(gene_mat))

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

# 5. 保存
saveRDS(gene_mat, "03_results/V172_Fixed/GSE59867_gene_level_FIXED.rds")
if(sum(found)>0) {
  extracted <- gene_mat[targets[found], , drop=FALSE]
  write.csv(data.frame(Gene=rownames(extracted), extracted, check.names=FALSE),
            "03_results/V172_Fixed/GSE59867_target_genes.csv", row.names=FALSE)
}

# 6. NDUFB7基础统计
if("NDUFB7" %in% rownames(gene_mat)) {
  vals <- gene_mat["NDUFB7",]
  message("[FOUND] NDUFB7: range ", round(min(vals),3), " - ", round(max(vals),3), 
          ", mean ", round(mean(vals),3))
}

message("[DONE] Saved to 03_results/V172_Fixed/")
