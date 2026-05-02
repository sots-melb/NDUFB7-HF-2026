# 铁死亡评分最终版
suppressPackageStartupMessages(library(data.table))

log_file <- Sys.getenv("R_LOG")
project_dir <- Sys.getenv("PROJECT")

cat("=== 开始 ===\n", file=log_file, append=TRUE)

expr_file <- file.path(project_dir, "03_results/02_tables/GSE57338_expression_matrix_annotated.csv")
expr <- fread(expr_file, check.names=FALSE)
cat("表达矩阵:", nrow(expr), "x", ncol(expr), "\n", file=log_file, append=TRUE)

# 核心铁死亡基因（与NDUFB7显著相关）
core_genes <- c("FTL", "SAT1", "NFE2L2", "SLC7A11", "ACSL4")
found <- intersect(core_genes, expr$gene_symbol)
cat("找到基因:", paste(found, collapse=", "), "\n", file=log_file, append=TRUE)

ferro_expr <- expr[gene_symbol %in% found]
ferro_mat <- as.matrix(ferro_expr[, -c("probe_id", "gene_symbol"), with=FALSE])
rownames(ferro_mat) <- ferro_expr$gene_symbol

# 标准化：行（基因）标准化，然后列（样本）均值
ferro_scaled <- t(scale(t(ferro_mat)))
ferro_score <- colMeans(ferro_scaled, na.rm=TRUE)

cat("评分样本数:", length(ferro_score), "\n", file=log_file, append=TRUE)
cat("评分范围:", round(range(ferro_score, na.rm=TRUE), 3), "\n", file=log_file, append=TRUE)

# NDUFB7相关性
ndufb7_file <- file.path(project_dir, "03_results/02_figures_tables/GSE57338_NDUFB7_expression_v3.csv")
ndufb7 <- fread(ndufb7_file, check.names=FALSE)
value_cols <- setdiff(names(ndufb7), c("probe_id", "gene_symbol", "ID_REF"))
ndufb7_vec <- as.numeric(ndufb7[1, ..value_cols])
names(ndufb7_vec) <- value_cols

common <- intersect(names(ferro_score), names(ndufb7_vec))
cat("共同样本:", length(common), "\n", file=log_file, append=TRUE)

ndufb7_match <- ndufb7_vec[common]
score_match <- ferro_score[common]

cor_test <- cor.test(ndufb7_match, score_match, method="pearson")
cat("Pearson r:", round(cor_test$estimate, 4), "p:", format(cor_test$p.value, scientific=TRUE, digits=2), "\n", file=log_file, append=TRUE)

# 保存
result <- data.table(n_samples=length(common), rho=cor_test$estimate, p_value=cor_test$p.value, genes=paste(found, collapse=","), timestamp=Sys.time())
fwrite(result, file.path(project_dir, "03_results/02_figures_tables/GSE57338_ferroptosis_FINAL.csv"))

sample_scores <- data.table(sample_id=names(ferro_score), ferroptosis_score=ferro_score)
fwrite(sample_scores, file.path(project_dir, "03_results/02_figures_tables/GSE57338_ferroptosis_sample_FINAL.csv"))

cat("=== 完成 ===\n", file=log_file, append=TRUE)
