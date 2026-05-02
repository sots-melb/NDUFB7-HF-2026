# --- 自动安装依赖 ---
for (pkg in c("lme4", "lmerTest", "UCell")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos="http://cran.us.r-project.org")
}
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(lme4))
suppressMessages(library(lmerTest))
suppressMessages(library(UCell))

message("▶ [1/4] 载入大对象并进行 20k 抽稀与心肌提纯...")
rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_figures_tables/GSE183852_scored.RDS"
if(!file.exists(rds_path)) stop("❌ 找不到 GSE183852_scored.RDS")

full_obj <- readRDS(rds_path)
set.seed(42)
keep_cells <- sample(Cells(full_obj), min(20000, ncol(full_obj)))
obj <- subset(full_obj, cells = keep_cells)
rm(full_obj); gc()

cm_cells <- WhichCells(obj, expression = TNNT2 > 0 & MYH7 > 0)
obj_cm <- subset(obj, cells = cm_cells)
message(sprintf("✅ 提纯完成，纯净 CM 数量: %d", ncol(obj_cm)))

# --- 安全打分函数 ---
safe_score <- function(seurat_obj, gene_list, score_name) {
    valid_genes <- intersect(gene_list, rownames(seurat_obj))
    if(length(valid_genes) > 0) {
        expr_matrix <- GetAssayData(seurat_obj, layer = "data")[valid_genes, , drop=FALSE]
        seurat_obj[[score_name]] <- colMeans(expr_matrix)
    } else {
        seurat_obj[[score_name]] <- 0
    }
    return(seurat_obj)
}

obj_cm <- safe_score(obj_cm, c('ATP5F1A', 'COX4I1', 'CYC1', 'NDUFA4', 'SDHA'), "OXPHOS_Score")
obj_cm <- safe_score(obj_cm, c('SLC7A11', 'GPX4', 'FTH1', 'FTL', 'NFE2L2', 'GSS'), "Ferro_Def_Score")

# --- 保存对象供以后画图用 ---
out_rds <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS"
saveRDS(obj_cm, out_rds)
message("✅ 对象已成功保存至硬盘，不再报错！")

# --- 提取数据用于稳健性分析 ---
ndufb7_expr <- as.numeric(GetAssayData(obj_cm, layer = "data")["NDUFB7", ])
df <- data.frame(
    NDUFB7 = ndufb7_expr,
    OXPHOS = obj_cm$OXPHOS_Score,
    Ferro_Def = obj_cm$Ferro_Def_Score,
    SampleID = if("orig.ident" %in% colnames(obj_cm@meta.data)) obj_cm$orig.ident else as.factor(obj_cm@meta.data[[1]])
)

message("\n=======================================================")
message("🛡️ 防御策略 1：Bootstrapping (100次重抽样 95% CI) 🛡️")
message("=======================================================")
set.seed(123)
boot_ox <- numeric(100); boot_fd <- numeric(100)
for(i in 1:100) {
    idx <- sample(1:nrow(df), replace = TRUE)
    boot_ox[i] <- cor(df$NDUFB7[idx], df$OXPHOS[idx], method="pearson")
    boot_fd[i] <- cor(df$NDUFB7[idx], df$Ferro_Def[idx], method="pearson")
}
ci_ox <- quantile(boot_ox, c(0.025, 0.975))
ci_fd <- quantile(boot_fd, c(0.025, 0.975))
message(sprintf("🌟 OXPHOS 产能:   r_mean = %.4f, 95%% CI [%.4f, %.4f]", mean(boot_ox), ci_ox[1], ci_ox[2]))
message(sprintf("🌟 铁死亡防线:    r_mean = %.4f, 95%% CI [%.4f, %.4f]", mean(boot_fd), ci_fd[1], ci_fd[2]))

message("\n=======================================================")
message("🛡️ 防御策略 2：线性混合效应模型 (LMM 排除批次效应) 🛡️")
message("=======================================================")
lmm_ox <- lmer(OXPHOS ~ NDUFB7 + (1|SampleID), data = df)
lmm_fd <- lmer(Ferro_Def ~ NDUFB7 + (1|SampleID), data = df)
summ_ox <- summary(lmm_ox)$coefficients
summ_fd <- summary(lmm_fd)$coefficients
message(sprintf("🌟 排除患者差异后，NDUFB7 对 OXPHOS 的固定效应: t-value = %.2f, p-value = %.2e", summ_ox["NDUFB7","t value"], summ_ox["NDUFB7","Pr(>|t|)"]))
message(sprintf("🌟 排除患者差异后，NDUFB7 对 铁死亡防线 的固定效应: t-value = %.2f, p-value = %.2e", summ_fd["NDUFB7","t value"], summ_fd["NDUFB7","Pr(>|t|)"]))

message("\n=======================================================")
message("🛡️ 防御策略 3：UCell 算法交叉验证 (非参数秩次打分) 🛡️")
message("=======================================================")
gene_list <- list(
  UCell_OXPHOS = intersect(c('ATP5F1A', 'COX4I1', 'CYC1', 'NDUFA4', 'SDHA'), rownames(obj_cm)),
  UCell_FerroDef = intersect(c('SLC7A11', 'GPX4', 'FTH1', 'FTL', 'NFE2L2', 'GSS'), rownames(obj_cm))
)
expr_mat <- GetAssayData(obj_cm, layer="data")
ucell_scores <- ScoreSignatures_UCell(expr_mat, features=gene_list, ncores=1)
r_ucell_ox <- cor.test(ndufb7_expr, ucell_scores[,"UCell_OXPHOS"])
r_ucell_fd <- cor.test(ndufb7_expr, ucell_scores[,"UCell_FerroDef"])
message(sprintf("🌟 UCell_OXPHOS 产能: r = %.4f (p=%.2e)", r_ucell_ox$estimate, r_ucell_ox$p.value))
message(sprintf("🌟 UCell_铁死亡防线:  r = %.4f (p=%.2e)", r_ucell_fd$estimate, r_ucell_fd$p.value))

message("\n🎉 完美通关！所有审稿人防御装甲已就绪！")
