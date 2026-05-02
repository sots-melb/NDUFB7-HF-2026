# 铁死亡评分修复脚本
suppressPackageStartupMessages({library(data.table)})

log_file <- Sys.getenv("R_LOG")
project_dir <- Sys.getenv("PROJECT")

cat("=== 铁死亡评分修复开始 ===\n", file=log_file, append=TRUE)

expr_file <- file.path(project_dir, "03_results/02_tables/GSE57338_expression_matrix_annotated.csv")
expr <- fread(expr_file, check.names=FALSE)
cat("表达矩阵维度:", nrow(expr), "x", ncol(expr), "\n", file=log_file, append=TRUE)

# 铁死亡基因集
ferro_genes <- c("GPX4", "SLC7A11", "FTH1", "FTL", "NCOA4", "ACSL4", "LPCAT3", "ALOX15", "SAT1", "NFE2L2", "SLC25A28", "KEAP1", "GSS", "GCLC")
found_genes <- intersect(ferro_genes, expr$gene_symbol)
cat("找到铁死亡基因:", length(found_genes), "/", length(ferro_genes), "\n", file=log_file, append=TRUE)
cat("具体基因:", paste(found_genes, collapse=", "), "\n", file=log_file, append=TRUE)

# 提取铁死亡基因表达矩阵
ferro_expr <- expr[gene_symbol %in% found_genes]
ferro_mat <- as.matrix(ferro_expr[, -c("probe_id", "gene_symbol"), with=FALSE])
rownames(ferro_mat) <- ferro_expr$gene_symbol
cat("铁死亡矩阵:", nrow(ferro_mat), "基因 x", ncol(ferro_mat), "样本\n", file=log_file, append=TRUE)

# 正向调控因子（促进铁死亡）
positive <- c("ACSL4", "LPCAT3", "ALOX15", "SAT1")
# 负向调控因子（抑制铁死亡）
negative <- c("GPX4", "SLC7A11", "FTH1", "FTL", "NFE2L2")

pos_in <- intersect(positive, found_genes)
neg_in <- intersect(negative, found_genes)
cat("正向基因:", paste(pos_in, collapse=", "), "\n", file=log_file, append=TRUE)
cat("负向基因:", paste(neg_in, collapse=", "), "\n", file=log_file, append=TRUE)

# === 核心修复: 行标准化（每个基因z-score），然后列均值（每个样本）===
sample_scores <- data.table(sample_id = colnames(ferro_mat))

if(length(pos_in) > 0) {
    pos_mat <- ferro_mat[pos_in, , drop=FALSE]
    # t(scale(t(mat))) = 对每行（基因）标准化，保持矩阵维度不变
    pos_scaled <- t(scale(t(pos_mat)))
    pos_scores <- colMeans(pos_scaled, na.rm=TRUE)
    sample_scores$positive_score <- pos_scores[sample_scores$sample_id]
    cat("正向评分样本数:", length(pos_scores), "范围:", round(range(pos_scores, na.rm=TRUE), 3), "\n", file=log_file, append=TRUE)
}

if(length(neg_in) > 0) {
    neg_mat <- ferro_mat[neg_in, , drop=FALSE]
    neg_scaled <- t(scale(t(neg_mat)))
    neg_scores <- colMeans(neg_scaled, na.rm=TRUE)
    sample_scores$negative_score <- neg_scores[sample_scores$sample_id]
    cat("负向评分样本数:", length(neg_scores), "范围:", round(range(neg_scores, na.rm=TRUE), 3), "\n", file=log_file, append=TRUE)
}

# 综合评分 = 正向 - 负向
if(all(c("positive_score", "negative_score") %in% names(sample_scores))) {
    sample_scores$composite_score <- sample_scores$positive_score - sample_scores$negative_score
} else if("positive_score" %in% names(sample_scores)) {
    sample_scores$composite_score <- sample_scores$positive_score
} else {
    sample_scores$composite_score <- -sample_scores$negative_score
}

cat("综合评分范围:", round(range(sample_scores$composite_score, na.rm=TRUE), 3), "\n", file=log_file, append=TRUE)

# === NDUFB7相关性 ===
ndufb7_file <- file.path(project_dir, "03_results/02_figures_tables/GSE57338_NDUFB7_expression_v3.csv")
if(file.exists(ndufb7_file)) {
    ndufb7 <- fread(ndufb7_file, check.names=FALSE)
    value_cols <- setdiff(names(ndufb7), c("probe_id", "gene_symbol", "ID_REF"))
    if(length(value_cols) > 0) {
        ndufb7_vec <- as.numeric(ndufb7[1, ..value_cols])
        names(ndufb7_vec) <- value_cols
        
        common_samples <- intersect(sample_scores$sample_id, names(ndufb7_vec))
        cat("共同样本数:", length(common_samples), "\n", file=log_file, append=TRUE)
        
        if(length(common_samples) > 10) {
            ndufb7_matched <- ndufb7_vec[common_samples]
            score_matched <- sample_scores[sample_id %in% common_samples, composite_score]
            
            cor_test <- cor.test(ndufb7_matched, score_matched, method="pearson")
            cat("NDUFB7-铁死亡评分 Pearson r:", round(cor_test$estimate, 4), 
                "p:", format(cor_test$p.value, scientific=TRUE, digits=2), "\n", file=log_file, append=TRUE)
            
            # 保存结果
            result <- data.table(
                n_samples = length(common_samples),
                rho = cor_test$estimate,
                p_value = cor_test$p.value,
                method = "pearson",
                positive_genes = paste(pos_in, collapse=","),
                negative_genes = paste(neg_in, collapse=","),
                timestamp = Sys.time()
            )
            fwrite(result, file.path(project_dir, "03_results/02_figures_tables/GSE57338_ferroptosis_scores_v2.csv"))
            cat("✅ 铁死亡评分结果已保存\n", file=log_file, append=TRUE)
            
            # 单基因相关性
            single_results <- list()
            for(g in found_genes) {
                if(g %in% rownames(ferro_mat)) {
                    g_vec <- ferro_mat[g, common_samples]
                    g_cor <- cor.test(ndufb7_matched, g_vec, method="pearson")
                    single_results[[g]] <- data.table(gene=g, rho=g_cor$estimate, p_value=g_cor$p.value)
                }
            }
            if(length(single_results) > 0) {
                single_df <- rbindlist(single_results)
                fwrite(single_df, file.path(project_dir, "03_results/02_figures_tables/GSE57338_ferroptosis_single_gene_cor.csv"))
                cat("✅ 单基因相关性已保存\n", file=log_file, append=TRUE)
            }
        }
    }
}

# 保存样本级评分
fwrite(sample_scores, file.path(project_dir, "03_results/02_figures_tables/GSE57338_ferroptosis_sample_scores.csv"))
cat("✅ 样本级评分已保存，行数:", nrow(sample_scores), "\n", file=log_file, append=TRUE)
cat("=== 铁死亡评分修复完成 ===\n", file=log_file, append=TRUE)
