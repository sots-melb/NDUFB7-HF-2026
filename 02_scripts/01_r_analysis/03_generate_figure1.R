#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     Figure 1 全套图表生成（A: 表达分布, B-D: 荟萃分析）  ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# 创建输出目录
dir.create("03_results/07_figures", showWarnings = FALSE, recursive = TRUE)

# ==================== Figure 1A: NDUFB7表达分布箱线图 ====================
cat("\n========== Figure 1A: 表达分布箱线图 ==========\n")

# 读取各数据集
gse57338 <- readRDS("01_data/01_raw_geo/GSE57338/GSE57338_gene_level.rds")$exprs
gse141910 <- readRDS("01_data/01_raw_geo/GSE141910/GSE141910_merged_matrix.rds")
gse5406 <- readRDS("01_data/01_raw_geo/GSE5406/GSE5406_gene_level.rds")$exprs
gse79962 <- readRDS("01_data/01_raw_geo/GSE79962/GSE79962_gene_level.rds")$exprs
gse116250 <- readRDS("01_data/01_raw_geo/GSE116250/GSE116250_rpkm_matrix.rds")$exprs

# 提取NDUFB7 + 分组信息
get_ndufb7_data <- function(gse_id, mat, ndufb7_id) {
  vals <- as.numeric(mat[ndufb7_id, ])
  samples <- colnames(mat)
  
  # 分组
  if (gse_id == "GSE116250") {
    grp <- ifelse(grepl("^NF", samples), "Control", "HF")
  } else if (gse_id == "GSE141910") {
    grp <- rep(NA, length(samples))
    grp[1:180] <- "HF"
    grp[181:366] <- "Control"
  } else {
    # 从已保存的分组文件读取
    groups_df <- read.csv("03_results/08_tables/TableS1_sample_groups_v3.csv")
    sub <- groups_df[groups_df$Dataset == gse_id, ]
    grp <- sub$Group[match(samples, sub$Sample)]
  }
  
  data.frame(
    Dataset = paste0(gse_id, "\n(", 
                    ifelse(gse_id=="GSE57338", "Affy ST, n=313",
                    ifelse(gse_id=="GSE141910", "RNA-seq, n=366",
                    ifelse(gse_id=="GSE5406", "U133A, n=210",
                    ifelse(gse_id=="GSE79962", "Affy ST, n=51",
                    "RNA-seq, n=64")))), ")"),
    Group = grp,
    Expression = vals,
    stringsAsFactors = FALSE
  )
}

plot_data <- rbind(
  get_ndufb7_data("GSE57338", gse57338, "NDUFB7"),
  get_ndufb7_data("GSE141910", gse141910, "ENSG00000167996"),
  get_ndufb7_data("GSE5406", gse5406, "NDUFB7"),
  get_ndufb7_data("GSE79962", gse79962, "NDUFB7"),
  get_ndufb7_data("GSE116250", gse116250, "ENSG00000167996")
)

# 绘制
pdf("03_results/07_figures/Fig1A_expression_boxplot.pdf", width = 12, height = 7)
par(mar = c(6, 4, 4, 2))
boxplot(Expression ~ Dataset + Group, data = plot_data, 
        main = "NDUFB7 Expression Across Datasets",
        ylab = "Expression Level (raw scale)",
        xlab = "",
        col = c("#FF6B6B", "#4ECDC4"),  # HF=红, Control=青
        las = 2, cex.axis = 0.7, cex.main = 1.2)
legend("topright", legend = c("HF", "Control"), 
       fill = c("#FF6B6B", "#4ECDC4"), title = "Group")
dev.off()
cat("[保存] Fig1A_expression_boxplot.pdf\n")

# ==================== Figure 1B: 主分析森林图（重新生成高清版） ====================
cat("\n========== Figure 1B: 主分析森林图 ==========\n")
library(metafor)

meta_data <- data.frame(
  id = c("GSE57338", "GSE141910", "GSE5406", "GSE79962", "GSE116250"),
  yi = c(0.07203215, -0.20456494, 0.30004097, -0.42079546, -0.98995972),
  vi = c(0.01301095, 0.01098907, 0.06786898, 0.11764506, 0.09908498),
  n = c(313, 366, 210, 51, 64),
  stringsAsFactors = FALSE
)

res <- rma(yi = yi, vi = vi, data = meta_data, method = "REML",
           slab = paste0(id, " (n=", n, ")"))

pdf("03_results/07_figures/Fig1B_forest_main.pdf", width = 11, height = 8)
forest(res, main = "NDUFB7 in Heart Failure: Pooled Analysis (k=5)",
       xlab = "Standardized Mean Difference (Cohen's d)",
       mlab = paste0("RE Model (I² = ", round(res$I2, 1), "%, p = ", format(res$QEp, digits = 3), ")"),
       cex = 0.9, col = "darkblue", border = "darkblue", addfit = TRUE, addpred = TRUE)
dev.off()
cat("[保存] Fig1B_forest_main.pdf\n")

# ==================== Figure 1E: 敏感性分析森林图（Leave-one-out） ====================
cat("\n========== Figure 1E: 敏感性分析图 ==========\n")

loo_data <- data.frame(
  study = c("Full model", "w/o GSE57338", "w/o GSE141910", "w/o GSE5406", 
            "w/o GSE79962", "w/o GSE116250"),
  d = c(-0.200, -0.298, -0.219, -0.311, -0.170, -0.049),
  lower = c(-0.585, -0.793, -0.764, -0.727, -0.641, -0.281),
  upper = c(0.184, 0.197, 0.326, 0.105, 0.302, 0.184),
  I2 = c(81.9, 77.4, 80.9, 82.8, 88.1, 49.6),
  stringsAsFactors = FALSE
)

pdf("03_results/07_figures/Fig1E_sensitivity.pdf", width = 10, height = 6)
par(mar = c(5, 8, 4, 2))
plot(0, 0, type = "n", xlim = c(-1, 0.5), ylim = c(1, nrow(loo_data)),
     xlab = "Pooled Cohen's d", ylab = "", yaxt = "n", main = "Sensitivity Analysis: Leave-One-Out")
axis(2, at = 1:nrow(loo_data), labels = loo_data$study, las = 1)
abline(v = 0, lty = 2, col = "gray")
for (i in 1:nrow(loo_data)) {
  color <- ifelse(i == 1, "darkblue", ifelse(i == nrow(loo_data), "red", "gray"))
  lwd <- ifelse(i %in% c(1, nrow(loo_data)), 2, 1)
  points(loo_data$d[i], i, pch = 18, col = color, cex = 1.5)
  lines(c(loo_data$lower[i], loo_data$upper[i]), c(i, i), col = color, lwd = lwd)
}
legend("bottomleft", legend = c("Full model", "w/o GSE116250", "Other"),
       col = c("darkblue", "red", "gray"), pch = 18, lwd = c(2, 2, 1))
dev.off()
cat("[保存] Fig1E_sensitivity.pdf\n")

cat("\n========================================\n")
cat("🎉 Figure 1 全套图表生成完成！\n")
cat("========================================\n")
cat("\n文件清单:\n")
files <- list.files("03_results/07_figures", pattern = "Fig1", full.names = TRUE)
for (f in files) {
  cat("  ", basename(f), "(", round(file.info(f)$size/1024, 1), "KB)\n")
}
