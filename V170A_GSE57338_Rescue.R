suppressPackageStartupMessages(library(data.table))
PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V170A: GSE57338 表达矩阵抢救")
message("========================================")

# 按优先级尝试所有可能的表达数据源
candidates <- list(
  list(file = "03_results/02_tables/GSE57338_expression_matrix_annotated.csv", type = "csv"),
  list(file = "03_results/02_tables/GSE57338_expression_matrix_raw.csv", type = "csv"),
  list(file = "03_results/03_pathway_analysis/hdWGCNA/gse57338_expr_matrix.rds", type = "rds"),
  list(file = "03_results/01_seurat_objects/GSE57338_gene_level.rds", type = "seurat"),
  list(file = "03_results/01_seurat_objects/GSE57338_pure.rds", type = "seurat")
)

expr_mat <- NULL
source_used <- NULL

for(cand in candidates) {
  f <- cand$file
  if(file.exists(f)) {
    message("\n[TRY] ", f, " (", cand$type, ")")
    tryCatch({
      if(cand$type == "csv") {
        # 尝试读取CSV，判断方向
        df <- fread(f, header = TRUE, nrows = 5)
        n_col <- ncol(df)
        n_row <- nrow(df)
        message("  Preview: ", n_row, " x ", n_col)
        message("  First col: ", names(df)[1])
        
        # 判断：第一列是基因名还是样本名
        first_vals <- as.character(df[[1]])
        looks_like_gene <- any(grepl("^NDUFB7$|^ENSG|^\\d+_at", first_vals, ignore.case = TRUE))
        
        if(looks_like_gene) {
          message("  [DIAG] Row = gene, Col = sample")
          # 重新读取全部
          df_full <- fread(f, header = TRUE)
          gene_col <- names(df_full)[1]
          samples <- setdiff(names(df_full), gene_col)
          expr_mat <- as.matrix(df_full[, ..samples])
          rownames(expr_mat) <- df_full[[gene_col]]
        } else {
          message("  [DIAG] Col = gene, Row = sample — transposing")
          df_full <- fread(f, header = TRUE)
          expr_mat <- as.matrix(t(df_full[, -1]))
          colnames(expr_mat) <- df_full[[1]]
        }
        source_used <- f
        break
      } else if(cand$type == "rds") {
        obj <- readRDS(f)
        if(is.matrix(obj)) {
          expr_mat <- obj
          source_used <- paste0(f, " (matrix)")
          message("  [PASS] Matrix: ", nrow(expr_mat), " x ", ncol(expr_mat))
        } else if(is.data.frame(obj)) {
          expr_mat <- as.matrix(obj)
          source_used <- paste0(f, " (data.frame)")
          message("  [PASS] Data.frame→Matrix: ", nrow(expr_mat), " x ", ncol(expr_mat))
        }
        break
      } else if(cand$type == "seurat") {
        obj <- readRDS(f)
        if(inherits(obj, "Seurat")) {
          suppressPackageStartupMessages(library(Seurat))
          expr_mat <- as.matrix(GetAssayData(obj, slot = "data"))
          source_used <- paste0(f, " (Seurat data slot)")
          message("  [PASS] Seurat: ", nrow(expr_mat), " x ", ncol(expr_mat))
        } else {
          message("  [SKIP] Not a Seurat object")
        }
        break
      }
    }, error = function(e) {
      message("  [FAIL] ", conditionMessage(e))
    })
  } else {
    message("\n[MISS] ", f)
  }
}

if(is.null(expr_mat)) {
  message("\n[FAIL] No valid expression matrix found in any candidate")
  # 保存诊断
  diag <- data.frame(
    Candidate = sapply(candidates, function(x) x$file),
    Exists = sapply(candidates, function(x) file.exists(x$file)),
    stringsAsFactors = FALSE
  )
  fwrite(diag, "03_results/V170_Rescue/V170A_diagnosis.csv")
  quit(save="no", status=1)
}

message("\n[PASS] Final matrix: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")
message("  Source: ", source_used)

# 保存标准化矩阵（基因×样本）
saveRDS(expr_mat, "03_results/V170_Rescue/GSE57338_rescued_matrix.rds")
fwrite(data.frame(Gene = rownames(expr_mat), expr_mat), 
       "03_results/V170_Rescue/GSE57338_rescued_matrix.csv")

# 验证NDUFB7
nduf_match <- grep("NDUFB7", rownames(expr_mat), value = TRUE, ignore.case = TRUE)
if(length(nduf_match) > 0) {
  message("\n[PASS] NDUFB7 found: ", paste(nduf_match, collapse = ", "))
  for(nm in nduf_match[1:min(3, length(nduf_match))]) {
    vals <- as.numeric(expr_mat[nm, ])
    message("  ", nm, " range: ", round(min(vals, na.rm = TRUE), 3), 
            " - ", round(max(vals, na.rm = TRUE), 3),
            " | median: ", round(median(vals, na.rm = TRUE), 3))
  }
} else {
  message("\n[WARN] NDUFB7 not in rownames")
  message("  Sample rownames: ", paste(head(rownames(expr_mat), 10), collapse = ", "))
}

message("\n[DONE] V170A")
