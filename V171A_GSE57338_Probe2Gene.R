suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V171A: GSE57338 探针→基因映射")
message("========================================")

# 1. 加载探针级矩阵（从V170A确认的位置）
rds_file <- "03_results/03_pathway_analysis/hdWGCNA/gse57338_expr_matrix.rds"
if(!file.exists(rds_file)) {
  message("[FAIL] RDS not found: ", rds_file)
  quit(status=1)
}

expr_mat <- readRDS(rds_file)
message("[PASS] Loaded matrix: ", nrow(expr_mat), " x ", ncol(expr_mat))
message("  Sample rownames: ", paste(head(rownames(expr_mat), 3), collapse=", "), "...")

# 2. 查找GPL11532注释文件
gpl_files <- c(
  "01_data/02_gpl_annotation/GPL11532.annot.gz",
  "01_data/02_gpl_annotation/GPL11532.annot",
  "Downloads/GPL11532.annot.gz",
  "Downloads/GPL11532.annot",
  "01_data/01_raw_geo/GPL11532/GPL11532.annot.gz"
)

gpl_file <- NULL
for(f in gpl_files) {
  if(file.exists(f)) {
    gpl_file <- f
    message("[FOUND] GPL annotation: ", f)
    break
  }
}

# 如果没有注释文件，尝试用series matrix中的注释
if(is.null(gpl_file)) {
  message("[WARN] No GPL11532 annotation file found. Trying series matrix...")
  sm_files <- list.files("01_data/01_raw_geo/GSE57338", pattern="series_matrix", full.names=TRUE)
  if(length(sm_files) > 0) {
    message("[TRY] Series matrix: ", sm_files[1])
    # 这里只是标记，实际解析在V171C中
  }
  message("[FAIL] Cannot map probes without annotation. Please download GPL11532.annot.gz")
  # 保存探针级矩阵供后续使用
  saveRDS(expr_mat, "03_results/V171_Fixed/GSE57338_probe_level.rds")
  write.csv(data.frame(ProbeID=rownames(expr_mat)), "03_results/V171_Fixed/GSE57338_probe_list.csv", row.names=FALSE)
  quit(status=1)
}

# 3. 解析GPL注释
message("[STEP 3] Parsing GPL annotation...")
gpl <- fread(gpl_file, skip="#ID", header=TRUE, fill=TRUE, quote="")

# 确定ID列和基因列
id_col <- grep("ID|id|Probe", names(gpl), value=TRUE)[1]
gene_col <- grep("Gene|Symbol|gene_symbol|GENE_SYMBOL", names(gpl), value=TRUE)[1]

if(is.na(id_col) || is.na(gene_col)) {
  message("[DIAG] Columns: ", paste(names(gpl)[1:5], collapse=", "))
  id_col <- names(gpl)[1]
  gene_col <- names(gpl)[which(grepl("Gene|Symbol", names(gpl), ignore.case=TRUE))[1]]
  if(is.na(gene_col)) gene_col <- names(gpl)[2]
}

message("  ID column: ", id_col)
message("  Gene column: ", gene_col)

# 创建映射表
mapping <- gpl[, c(id_col, gene_col), with=FALSE]
setnames(mapping, c("ProbeID", "GeneSymbol"))
mapping <- mapping[!is.na(GeneSymbol) & GeneSymbol != "" & GeneSymbol != "---"]
mapping <- mapping[!duplicated(ProbeID)]

message("[PASS] Mapping table: ", nrow(mapping), " probes")

# 4. 映射探针到基因
probe_ids <- rownames(expr_mat)
matched <- probe_ids %in% mapping$ProbeID
message("[STEP 4] Matching probes...")
message("  Matched: ", sum(matched), " / ", length(probe_ids))

# 提取可映射的探针
expr_mapped <- expr_mat[matched, , drop=FALSE]
gene_symbols <- mapping$GeneSymbol[match(rownames(expr_mapped), mapping$ProbeID)]

# 处理重复基因：取中位数
gene_df <- data.frame(GeneSymbol=gene_symbols, expr_mapped, check.names=FALSE)
gene_unique <- gene_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(everything(), median), .groups="drop")

# 转换为矩阵
gene_mat <- as.matrix(gene_unique[,-1])
rownames(gene_mat) <- gene_unique$GeneSymbol

message("[PASS] Gene-level matrix: ", nrow(gene_mat), " genes x ", ncol(gene_mat), " samples")

# 5. 检查NDUFB7
if("NDUFB7" %in% rownames(gene_mat)) {
  ndufb7_vals <- gene_mat["NDUFB7", ]
  message("[FOUND] NDUFB7 in gene matrix!")
  message("  Range: ", round(min(ndufb7_vals), 3), " - ", round(max(ndufb7_vals), 3))
  message("  Mean: ", round(mean(ndufb7_vals), 3))
} else {
  message("[MISS] NDUFB7 not in gene matrix")
  # 检查是否有近似名称
  ndufb7_like <- grep("NDUFB|NDUF", rownames(gene_mat), value=TRUE)
  if(length(ndufb7_like) > 0) {
    message("  Similar names: ", paste(head(ndufb7_like, 5), collapse=", "))
  }
}

# 6. 保存
saveRDS(gene_mat, "03_results/V171_Fixed/GSE57338_gene_level_FIXED.rds")
write.csv(data.frame(Gene=rownames(gene_mat), gene_mat, check.names=FALSE), 
          "03_results/V171_Fixed/GSE57338_gene_level_FIXED.csv", row.names=FALSE)

message("[DONE] Saved to 03_results/V171_Fixed/")
