#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(mclust)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V125_Distribution_Contrast")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V125_FINAL: 分布模式仲裁（暴力提取版）")
message("========================================")

# --- GSE183852: 暴力提取NDUFB7 ---
cm_file <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"
srt <- readRDS(cm_file)

# 暴力探测所有可能的counts位置
expr_183852 <- NULL
try_list <- list(
  tryCatch(as.numeric(srt@assays$RNA@counts["NDUFB7", ]), error=function(e) NULL),
  tryCatch(as.numeric(srt@assays$RNA$counts["NDUFB7", ]), error=function(e) NULL),
  tryCatch(as.numeric(srt@assays$RNA@layers$counts["NDUFB7", ]), error=function(e) NULL),
  tryCatch(as.numeric(srt@assays$RNA@data["NDUFB7", ]), error=function(e) NULL),
  tryCatch(as.numeric(srt[["RNA"]]$counts["NDUFB7", ]), error=function(e) NULL)
)
for (i in seq_along(try_list)) {
  if (!is.null(try_list[[i]]) && length(try_list[[i]]) > 100) {
    expr_183852 <- try_list[[i]]
    message("[FOUND] Method ", i, ": n=", length(expr_183852))
    break
  }
}
if (is.null(expr_183852)) {
  # 终极方案：遍历所有 assay slots
  for (slot in c("counts", "data", "scale.data")) {
    for (assay_name in names(srt@assays)) {
      mat <- tryCatch(slot(srt@assays[[assay_name]], slot), error=function(e) NULL)
      if (!is.null(mat) && "NDUFB7" %in% rownames(mat)) {
        expr_183852 <- as.numeric(mat["NDUFB7", ])
        message("[FOUND] Ultimate: assay=", assay_name, " slot=", slot, " n=", length(expr_183852))
        break
      }
    }
    if (!is.null(expr_183852)) break
  }
}
if (is.null(expr_183852)) stop("[FAIL] Cannot extract NDUFB7 from any slot")

expr_183852 <- expr_183852[!is.na(expr_183852)]
message("[PASS] GSE183852: n=", length(expr_183852), " zero_pct=", round(mean(expr_183852==0)*100,1), "%")

# --- GSE121893 ---
f121893 <- "~/Downloads/GSE121893_human_heart_sc_umi.csv.gz"
if (!file.exists(f121893)) stop("GSE121893 not found")
dt <- data.table::fread(f121893, header = TRUE)
gene_col <- colnames(dt)[1]
idx <- which(toupper(dt[[gene_col]]) == "NDUFB7")[1]
expr_121893 <- as.numeric(dt[idx, -1, with = FALSE])
expr_121893 <- expr_121893[!is.na(expr_121893)]
message("[PASS] GSE121893: n=", length(expr_121893), " zero_pct=", round(mean(expr_121893==0)*100,1), "%")

# 标准化
norm_183852 <- (expr_183852 - min(expr_183852)) / (max(expr_183852) - min(expr_183852) + 1e-10)
norm_121893 <- (expr_121893 - min(expr_121893)) / (max(expr_121893) - min(expr_121893) + 1e-10)

# BIC
bic_table <- data.frame()
for (g in 1:3) {
  mod_183852 <- tryCatch(densityMclust(norm_183852, G = g, verbose = FALSE), error = function(e) NULL)
  mod_121893 <- tryCatch(densityMclust(norm_121893, G = g, verbose = FALSE), error = function(e) NULL)
  bic_183852 <- if (!is.null(mod_183852)) mod_183852$bic else NA
  bic_121893 <- if (!is.null(mod_121893)) mod_121893$bic else NA
  bic_table <- rbind(bic_table, data.frame(Dataset=c("GSE183852","GSE121893"), G=g, BIC=c(bic_183852,bic_121893)))
}
best_183852 <- bic_table$G[which.max(bic_table$BIC[bic_table$Dataset=="GSE183852"])]
best_121893 <- bic_table$G[which.max(bic_table$BIC[bic_table$Dataset=="GSE121893"])]

write.csv(bic_table, file.path(outdir, "V125_BIC_model_selection.csv"), row.names = FALSE)
message("\nBest G: GSE183852=", best_183852, " GSE121893=", best_121893)

# KS + Fisher
ks_test <- ks.test(norm_183852, norm_121893)
fisher_mat <- matrix(c(sum(expr_183852==0), sum(expr_183852>0), sum(expr_121893==0), sum(expr_121893>0)), nrow=2)
fisher_test <- fisher.test(fisher_mat)

# 密度图
df_plot <- data.frame(Expression=c(norm_183852, norm_121893), 
                      Dataset=rep(c("GSE183852 (n=637)","GSE121893 (n=4933)"), c(length(norm_183852), length(norm_121893))))
p <- ggplot(df_plot, aes(x=Expression, fill=Dataset)) + geom_density(alpha=0.4, linewidth=0.8) +
  scale_fill_manual(values=c("GSE183852 (n=637)"="#440154","GSE121893 (n=4933)"="#35B779")) +
  annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5, size=3,
           label=paste0("KS D=",round(ks_test$statistic,3),"\np=",format(ks_test$p.value,digits=1,scientific=TRUE),
                       "\nBest G=",best_183852," vs ",best_121893)) +
  labs(title="NDUFB7 Distribution: GSE183852 vs GSE121893", x="Normalized Expression", y="Density") + theme_minimal()
ggsave(file.path(outdir,"V125_density_comparison.png"), p, width=7, height=5, dpi=300)

message("\n[DONE] V125_FINAL: ", outdir)
