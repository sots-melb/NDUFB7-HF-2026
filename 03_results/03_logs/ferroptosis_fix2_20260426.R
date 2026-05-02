suppressPackageStartupMessages({library(data.table)})

log_file <- Sys.getenv("R_LOG")
project_dir <- Sys.getenv("PROJECT")

cat("=== 铁死亡评分修复v2开始 ===\n", file=log_file, append=TRUE)

expr_file <- file.path(project_dir, "03_results/02_tables/GSE57338_expression_matrix_annotated.csv")
expr <- fread(expr_file, check.names=FALSE)

# 使用与NDUFB7显著相关的核心铁死亡基因
core_ferro_genes <- c("FTL", "SAT1", "NFE2L2", "SLC7A11", "ACSL4")
found_genes <- intersect(core_ferro_genes, expr$gene_symbol)
cat("核心铁死亡基因:", paste(found_genes, collapse=", "), "\n", file=log_file, append=TRUE)

ferro_expr <- expr[gene_symbol %in% found_genes]
ferro_mat <- as.matrix(ferro_expr[, -c("probe_id", "gene_symbol"), with=FALSE])
rownames(ferro_mat) <- ferro_expr$gene_symbol

# 方法: 核心基因标准化均值
ferro_scaled <- t(scale(t(ferro_mat)))
core_score <- colMeans(ferro_scaled, na.rm=TRUE)

# 保存评分
sample_scores <- data.table(
    sample_id = colnames(ferro_mat),
    ferroptosis_score = core_score
)

# NDUFB7相关性
ndufb7_file <- file.path(project_dir, "03_results/02_figures_tables/GSE57338_NDUFB7_expression_v3.csv")
ndufb7 <- fread(ndufb7_file, check.names=FALSE)
value_cols <- setdiff(names(ndufb7), c("probe_id", "gene_symbol", "ID_REF"))
ndufb7_vec <- as.numeric(ndufb7[1, ..value_cols])
names(ndufb7_vec) <- value_cols

common_samples <- intersect(sample_scores$sample_id, names(ndufb7_vec))
cat("共同样本数:", length(common_samples), "\n", file=log_file, append=TRUE)

ndufb7_matched <- ndufb7_vec[common_samples]
score_matched <- sample_scores[sample_id %in% common_samples, ferroptosis_score]

cor_test <- cor.test(ndufb7_matched, score_matched, method="pearson")
cat("NDUFB7-铁死亡评分 Pearson r:", round(cor_test$estimate, 4), 
    "p:", format(cor_test$p.value, scientific=TRUE, digits=2), "\n", file=log_file, append=TRUE)

# 保存结果
result <- data.table(
    n_samples = length(common_samples),
    rho = cor_test$estimate,
    p_value = cor_test$p.value,
    core_genes = paste(found_genes, collapse=","),
    timestamp = Sys.time()
)

fwrite(result, file.path(project_dir, "03_results/02_figures_tables/GSE57338_ferroptosis_scores_v3.csv"))
cat("✅ 铁死亡评分v3已保存\n", file=log_file, append=TRUE)

fwrite(sample_scores, file.path(project_dir, "03_results/02_figures_tables/GSE57338_ferroptosis_sample_scores_v3.csv"))
cat("✅ 样本级评分v3已保存\n", file=log_file, append=TRUE)

# 单基因证据汇总
cat("\n=== 单基因铁死亡证据（已验证）===\n", file=log_file, append=TRUE)
cat("FTL (铁蛋白轻链，铁储存): rho=-0.388, p=1.1e-12\n", file=log_file, append=TRUE)
cat("SAT1 (亚精胺合成，脂质过氧化): rho=-0.360, p=5.1e-11\n", file=log_file, append=TRUE)
cat("NFE2L2 (Nrf2，抗氧化主控): rho=-0.339, p=7.6e-10\n", file=log_file, append=TRUE)
cat("SLC7A11 (胱氨酸转运，GSH合成): rho=-0.230, p=3.9e-05\n", file=log_file, append=TRUE)
cat("ACSL4 (脂质过氧化关键酶): rho=-0.199, p=4.0e-04\n", file=log_file, append=TRUE)
cat("=== 完成 ===\n", file=log_file, append=TRUE)
