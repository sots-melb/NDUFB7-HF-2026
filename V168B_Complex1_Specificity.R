#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V168_PanDeath_Solidity/V168B_Complex1")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V168B: Complex I特异性验证")
message("========================================")

# === 1. 读取GSE57338表达矩阵 ===
expr_file <- "03_results/02_tables/GSE57338_expression_matrix_annotated.csv"
if(!file.exists(expr_file)){
  # 回退到raw
  expr_file <- "03_results/02_tables/GSE57338_expression_matrix_raw.csv"
}

if(!file.exists(expr_file)){
  message("[FAIL] GSE57338 expression matrix not found")
  quit(save="no", status=1)
}

message("\n[1/4] 读取表达矩阵: ", basename(expr_file))
df <- fread(expr_file, header = TRUE)
message("  Dimensions: ", nrow(df), " × ", ncol(df))
message("  第一列: ", names(df)[1])

# 判断方向：行=基因还是列=基因
gene_col <- names(df)[1]
is_gene_row <- any(grepl("NDUFB7|ENSG|Gene", df[[1]][1:5], ignore.case = TRUE))

if(is_gene_row){
  message("  [DIAG] 行=基因，列=样本")
  # 提取Complex I基因
  c1_genes <- c("NDUFB7", "NDUFA13", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB8", 
                "NDUFB9", "NDUFB10", "NDUFB11", "NDUFA1", "NDUFA2", "NDUFA3",
                "NDUFS1", "NDUFS2", "NDUFS3", "NDUFV1", "NDUFV2")
  
  # 找基因行
  find_gene <- function(g){
    idx <- which(tolower(df[[gene_col]]) == tolower(g))
    if(length(idx) == 0) idx <- grep(paste0("^", g, "$"), df[[gene_col]], ignore.case = TRUE)
    if(length(idx) == 0) idx <- grep(g, df[[gene_col]], ignore.case = TRUE)[1]
    return(idx)
  }
  
  c1_idx <- sapply(c1_genes, find_gene)
  names(c1_idx) <- c1_genes
  
  found_genes <- c1_genes[!is.na(c1_idx)]
  missing_genes <- c1_genes[is.na(c1_idx)]
  message("  Complex I基因找到: ", length(found_genes), "/", length(c1_genes))
  if(length(missing_genes) > 0) message("  缺失: ", paste(missing_genes, collapse = ", "))
  
  # 提取表达值
  num_cols <- setdiff(names(df), gene_col)
  expr_mat <- as.matrix(df[unlist(c1_idx[!is.na(c1_idx)]), ..num_cols])
  rownames(expr_mat) <- found_genes
  expr_mat <- apply(expr_mat, 2, as.numeric)
  
} else {
  message("  [DIAG] 列=基因，需要转置")
  # 寻找NDUFB7列
  c1_genes <- c("NDUFB7", "NDUFA13", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB8", 
                "NDUFB9", "NDUFB10", "NDUFB11", "NDUFA1", "NDUFA2", "NDUFA3")
  c1_cols <- sapply(c1_genes, function(g) grep(paste0("^", g, "$"), names(df), ignore.case = TRUE)[1])
  found_genes <- c1_genes[!is.na(c1_cols)]
  expr_mat <- t(as.matrix(df[, unlist(c1_cols[!is.na(c1_cols)]), with = FALSE]))
  rownames(expr_mat) <- found_genes
}

message("  表达矩阵: ", nrow(expr_mat), " genes × ", ncol(expr_mat), " samples")

# === 2. 读取或计算死亡通路评分 ===
score_candidates <- c(
  "03_results/V83_pyroptosis_main/V83B_GSE57338_sample_scores.csv",
  "03_results/02_figures_tables/GSE57338_ferroptosis_sample_scores_v3.csv",
  "03_results/02_figures_tables/GSE57338_ferroptosis_sample_scores.csv"
)

scores <- NULL
for(f in score_candidates){
  if(file.exists(f)){
    message("\n[2/4] 读取样本评分: ", basename(f))
    scores <- fread(f)
    message("  Columns: ", paste(names(scores)[1:min(10, ncol(scores))], collapse = ", "))
    if(ncol(scores) > 10) message("  ... and ", ncol(scores)-10, " more")
    break
  }
}

# 如果没有样本评分文件，从表达矩阵直接计算
if(is.null(scores)){
  message("\n[2/4] 无样本评分文件，从表达矩阵直接计算通路评分...")
  # 定义通路基因
  pw_genes <- list(
    Apoptosis = c("BCL2", "BAX", "CASP3", "CASP7", "CASP9", "BAD", "BAK1"),
    Necroptosis = c("RIPK1", "RIPK3", "MLKL", "CYLD", "TNF"),
    Ferroptosis_Defense = c("GPX4", "SLC7A11", "FTH1", "FTL", "NFE2L2"),
    Ferroptosis_Execution = c("ACSL4", "LPCAT3", "LOX", "PTGS2"),
    Autophagy = c("ATG5", "ATG7", "BECN1", "MAP1LC3B", "SQSTM1"),
    Pyroptosis = c("NLRP3", "GSDMD", "IL1B", "IL18", "CASP1")
  )
  
  # 在表达矩阵中找这些基因
  all_genes <- df[[gene_col]]
  sample_names <- setdiff(names(df), gene_col)
  
  scores_list <- list()
  for(pw in names(pw_genes)){
    g_found <- pw_genes[[pw]][tolower(pw_genes[[pw]]) %in% tolower(all_genes)]
    if(length(g_found) > 0){
      idx <- match(tolower(g_found), tolower(all_genes))
      vals <- as.matrix(df[idx, ..sample_names])
      vals <- apply(vals, 2, as.numeric)
      scores_list[[pw]] <- colMeans(vals, na.rm = TRUE)
      message("  ", pw, ": ", length(g_found), "/", length(pw_genes[[pw]], " genes found"))
    }
  }
  
  if(length(scores_list) > 0){
    scores <- as.data.table(scores_list)
    scores$Sample <- sample_names
  } else {
    message("[FAIL] No pathway genes found in expression matrix")
    quit(save="no", status=1)
  }
}

# === 3. 计算Complex I基因与通路评分的相关 ===
message("\n[3/4] 计算Complex I基因与死亡通路的相关...")

# 确保样本名匹配
common_samples <- intersect(colnames(expr_mat), names(scores))
if(length(common_samples) < 50){
  # 尝试样本名转换
  message("  [WARN] 样本名不匹配，尝试模糊匹配...")
  # scores的样本名可能是行名或列
  if("Sample" %in% names(scores)){
    score_samples <- scores$Sample
  } else {
    score_samples <- names(scores)
  }
  # 简化匹配逻辑：如果expr_mat列名含GSM，scores也应含GSM
  common_samples <- intersect(colnames(expr_mat), score_samples)
}

if(length(common_samples) < 50){
  message("[FAIL] Only ", length(common_samples), " common samples")
  # 保存诊断信息
  fwrite(data.frame(Expr_Samples = head(colnames(expr_mat), 5), Score_Samples = head(score_samples, 5)), 
         file.path(outdir, "V168B_sample_name_mismatch.csv"))
  quit(save="no", status=1)
}

message("  Common samples: ", length(common_samples))

expr_sub <- expr_mat[, common_samples]
scores_sub <- scores[match(common_samples, scores$Sample), ]

# 如果scores是data.table且Sample是列
if("Sample" %in% names(scores)){
  scores_sub <- scores[Sample %in% common_samples]
  scores_mat <- as.matrix(scores_sub[, setdiff(names(scores_sub), "Sample"), with = FALSE])
} else {
  scores_mat <- as.matrix(scores_sub)
}

# 确保数值
scores_mat <- apply(scores_mat, 2, as.numeric)

# 计算Spearman相关
results <- data.frame()
for(g in rownames(expr_sub)){
  for(j in 1:ncol(scores_mat)){
    pw <- colnames(scores_mat)[j]
    r <- cor(expr_sub[g, ], scores_mat[,j], method = "spearman", use = "pairwise.complete.obs")
    if(is.finite(r)){
      results <- rbind(results, data.frame(
        Gene = g,
        Is_NDUFB7 = g == "NDUFB7",
        Pathway = pw,
        Spearman_Rho = round(r, 3),
        N = sum(is.finite(expr_sub[g, ]) & is.finite(scores_mat[,j])),
        stringsAsFactors = FALSE
      ))
    }
  }
}

fwrite(results, file.path(outdir, "V168B_complex1_death_correlation.csv"))
message("\n=== Complex I vs 死亡通路相关 ===")
print(results[order(abs(Spearman_Rho), decreasing = TRUE), ])

# === 4. NDUFB7特异性判定 ===
message("\n[4/4] NDUFB7特异性判定...")

ndufb7_res <- results[results$Gene == "NDUFB7", ]
others_res <- results[results$Gene != "NDUFB7", ]

summary <- data.frame()
for(pw in unique(results$Pathway)){
  nduf_rho <- ndufb7_res$Spearman_Rho[ndufb7_res$Pathway == pw]
  other_rho <- others_res$Spearman_Rho[others_res$Pathway == pw]
  
  if(length(nduf_rho) > 0 && length(other_rho) > 0){
    is_max <- abs(nduf_rho) >= max(abs(other_rho), na.rm = TRUE)
    summary <- rbind(summary, data.frame(
      Pathway = pw,
      NDUFB7_Rho = round(nduf_rho, 3),
      Max_Other_Rho = round(max(abs(other_rho), na.rm = TRUE), 3),
      Median_Other_Rho = round(median(abs(other_rho), na.rm = TRUE), 3),
      NDUFB7_Is_Max = is_max,
      stringsAsFactors = FALSE
    ))
  }
}

fwrite(summary, file.path(outdir, "V168B_ndufb7_specificity_summary.csv"))
message("\n=== NDUFB7特异性汇总 ===")
print(summary)

# 判定
n_max <- sum(summary$NDUFB7_Is_Max, na.rm = TRUE)
message("\n[RESULT] NDUFB7 is the STRONGEST correlate in ", n_max, "/", nrow(summary), " pathways")
if(n_max >= 4){
  message("[PASS] NDUFB7 shows SPECIFIC pan-death signature vs other Complex I genes")
} else if(n_max >= 2){
  message("[PARTIAL] NDUFB7 is among the top correlates")
} else {
  message("[WARN] Other Complex I genes have stronger correlations — NDUFB7 may not be specific")
}

message("\n[DONE] V168B: ", outdir)
