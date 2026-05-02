#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(mclust)
  library(dplyr)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V125_Distribution_Contrast")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V125: GSE121893 vs GSE183852 分布模式统计仲裁")
message("========================================")

# --- 加载GSE183852 CM ---
cm_file <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"
srt <- readRDS(cm_file)
expr_183852 <- as.numeric(FetchData(srt, vars = "NDUFB7")$NDUFB7)
expr_183852 <- expr_183852[!is.na(expr_183852)]
message("[PASS] GSE183852: n=", length(expr_183852))

# --- 加载GSE121893 ---
f121893 <- "~/Downloads/GSE121893_human_heart_sc_umi.csv.gz"
if (!file.exists(f121893)) stop("GSE121893 not found")
dt <- data.table::fread(f121893, header = TRUE)
gene_col <- colnames(dt)[1]
idx <- which(toupper(dt[[gene_col]]) == "NDUFB7")[1]
expr_121893 <- as.numeric(dt[idx, -1, with = FALSE])
expr_121893 <- expr_121893[!is.na(expr_121893)]
message("[PASS] GSE121893: n=", length(expr_121893))

# --- 统一尺度（标准化到0-1）---
norm_183852 <- (expr_183852 - min(expr_183852)) / (max(expr_183852) - min(expr_183852) + 1e-10)
norm_121893 <- (expr_121893 - min(expr_121893)) / (max(expr_121893) - min(expr_121893) + 1e-10)

# --- BIC模型选择 ---
message("\n>>> [1/4] BIC模型选择（G=1,2,3）")

bic_table <- data.frame()
for (g in 1:3) {
  mod_183852 <- tryCatch(densityMclust(norm_183852, G = g, verbose = FALSE), error = function(e) NULL)
  mod_121893 <- tryCatch(densityMclust(norm_121893, G = g, verbose = FALSE), error = function(e) NULL)
  
  bic_183852 <- if (!is.null(mod_183852)) mod_183852$bic else NA
  bic_121893 <- if (!is.null(mod_121893)) mod_121893$bic else NA
  
  bic_table <- rbind(bic_table, data.frame(
    Dataset = c("GSE183852", "GSE121893"),
    G = g,
    BIC = c(bic_183852, bic_121893),
    Winner = c(ifelse(g == which.max(c(NA,NA,NA)), "?", ""), "")  # 稍后填充
  ))
}

# 确定最佳G
best_183852 <- bic_table$G[which.max(bic_table$BIC[bic_table$Dataset == "GSE183852"])]
best_121893 <- bic_table$G[which.max(bic_table$BIC[bic_table$Dataset == "GSE121893"])]
bic_table$Winner[bic_table$Dataset == "GSE183852" & bic_table$G == best_183852] <- "BEST"
bic_table$Winner[bic_table$Dataset == "GSE121893" & bic_table$G == best_121893] <- "BEST"

write.csv(bic_table, file.path(outdir, "V125_BIC_model_selection.csv"), row.names = FALSE)
message("GSE183852 best G: ", best_183852, " (BIC=", round(max(bic_table$BIC[bic_table$Dataset=="GSE183852"], na.rm=TRUE), 1), ")")
message("GSE121893 best G: ", best_121893, " (BIC=", round(max(bic_table$BIC[bic_table$Dataset=="GSE121893"], na.rm=TRUE), 1), ")")

# --- Kolmogorov-Smirnov检验 ---
message("\n>>> [2/4] 分布差异检验")
ks_test <- ks.test(norm_183852, norm_121893)
message("KS test D = ", round(ks_test$statistic, 4), ", p = ", format(ks_test$p.value, digits = 2, scientific = TRUE))

# --- 零值比例检验（Fisher精确检验）---
zero_183852 <- sum(expr_183852 == 0)
nonzero_183852 <- sum(expr_183852 > 0)
zero_121893 <- sum(expr_121893 == 0)
nonzero_121893 <- sum(expr_121893 > 0)

fisher_mat <- matrix(c(zero_183852, nonzero_183852, zero_121893, nonzero_121893), nrow = 2)
fisher_test <- fisher.test(fisher_mat)
message("\nZero-inflation Fisher exact test: OR=", round(fisher_test$estimate, 2), 
        ", p=", format(fisher_test$p.value, digits = 2, scientific = TRUE))

# --- Bootstrap BIC稳定性 ---
message("\n>>> [3/4] Bootstrap BIC稳定性（80%子抽样，100次）")

boot_bic <- function(expr, n_boot = 100) {
  sapply(1:n_boot, function(i) {
    samp <- sample(expr, size = floor(0.8 * length(expr)), replace = FALSE)
    sapply(1:3, function(g) {
      mod <- tryCatch(densityMclust(samp, G = g, verbose = FALSE), error = function(e) NULL)
      if (!is.null(mod)) mod$bic else NA
    })
  })
}

message("  Bootstrapping GSE183852...")
boot_183852 <- boot_bic(norm_183852)
best_g_boot_183852 <- apply(boot_183852, 2, function(x) which.max(x))
message("  GSE183852 best-G frequency: G1=", sum(best_g_boot_183852==1), 
        " G2=", sum(best_g_boot_183852==2), " G3=", sum(best_g_boot_183852==3))

message("  Bootstrapping GSE121893...")
boot_121893 <- boot_bic(norm_121893)
best_g_boot_121893 <- apply(boot_121893, 2, function(x) which.max(x))
message("  GSE121893 best-G frequency: G1=", sum(best_g_boot_121893==1), 
        " G2=", sum(best_g_boot_121893==2), " G3=", sum(best_g_boot_121893==3))

# --- 可视化：密度叠加 ---
message("\n>>> [4/4] 密度对比图")
df_plot <- data.frame(
  Expression = c(norm_183852, norm_121893),
  Dataset = rep(c("GSE183852 (n=637)", "GSE121893 (n=4933)"), c(length(norm_183852), length(norm_121893)))
)

p <- ggplot(df_plot, aes(x = Expression, fill = Dataset)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  scale_fill_manual(values = c("GSE183852 (n=637)" = "#440154", "GSE121893 (n=4933)" = "#35B779")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3,
           label = paste0("KS D=", round(ks_test$statistic, 3), 
                         "\np=", format(ks_test$p.value, digits=1, scientific=TRUE),
                         "\nGSE183852 best G=", best_183852,
                         "\nGSE121893 best G=", best_121893)) +
  labs(title = "NDUFB7 Distribution: GSE183852 vs GSE121893",
       subtitle = "Normalized to 0-1 scale | G=3 (tri-modal) vs G=2 (bi-modal)",
       x = "Normalized NDUFB7 Expression", y = "Density") +
  theme_minimal(base_size = 10)

ggsave(file.path(outdir, "V125_density_comparison.png"), p, width = 7, height = 5, dpi = 300)

# --- 保存统计汇总 ---
stats_summary <- data.frame(
  Test = c("KS_D", "KS_p", "Fisher_OR", "Fisher_p", "BIC_best_G_183852", "BIC_best_G_121893",
           "Boot_G3_freq_183852", "Boot_G2_freq_121893"),
  Value = c(ks_test$statistic, ks_test$p.value, fisher_test$estimate, fisher_test$p.value,
            best_183852, best_121893,
            sum(best_g_boot_183852==3)/length(best_g_boot_183852),
            sum(best_g_boot_121893==2)/length(best_g_boot_121893))
)
write.csv(stats_summary, file.path(outdir, "V125_statistical_summary.csv"), row.names = FALSE)

message("\n=== V125 统计汇总 ===")
print(stats_summary)

# --- 判定 ---
message("\n=== 判定 ===")
if (best_183852 == 3 && best_121893 == 2 && ks_test$p.value < 0.05) {
  message("[PASS] 分布模式差异有统计支持：GSE183852=G=3, GSE121893=G=2, 分布显著不同")
  message("[NARRATIVE] 'Platform- or etiology-dependent distribution modality'")
} else if (ks_test$p.value > 0.05) {
  message("[WARN] 分布差异不显著，可能为同一模式的采样变异性")
  message("[NARRATIVE] 降级为'样本量依赖性检测灵敏度差异'")
} else {
  message("[PARTIAL] BIC选择不一致，需谨慎解读")
}

message("[DONE] V125: ", outdir)
