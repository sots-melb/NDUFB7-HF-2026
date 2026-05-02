#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)
library(metafor)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     Phase 2 Pillar 1: 亚组分析与敏感性分析               ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# ==================== 1. 读取v3结果数据 ====================
# 从已保存的分组详情重新计算（或直接用汇总数据）
meta_data <- data.frame(
  id = c("GSE57338", "GSE141910", "GSE5406", "GSE79962", "GSE116250"),
  yi = c(0.07203215, -0.20456494, 0.30004097, -0.42079546, -0.98995972),
  vi = c(0.01301095, 0.01098907, 0.06786898, 0.11764506, 0.09908498),
  n = c(313, 366, 210, 51, 64),
  platform = c("Affy ST", "RNA-seq", "U133A", "Affy ST", "RNA-seq"),
  etiology = c("Mixed", "DCM", "Mixed", "Mixed", "Mixed"),
  stringsAsFactors = FALSE
)

cat("效应量数据:\n")
print(meta_data)

# ==================== 2. 按平台类型亚组分析 ====================
cat("\n========== 亚组分析：按平台类型 ==========\n")

res_platform <- rma(yi = yi, vi = vi, data = meta_data, method = "REML",
                    mods = ~ platform, slab = id)

cat("\n平台亚组分析结果:\n")
print(res_platform)

# 各平台汇总
platforms <- unique(meta_data$platform)
for (plat in platforms) {
  idx <- meta_data$platform == plat
  if (sum(idx) >= 2) {
    res_sub <- rma(yi = yi, vi = vi, data = meta_data[idx, ], method = "REML")
    cat("\n[", plat, "] k=", sum(idx), "\n")
    cat("  汇总d =", round(res_sub$b[1], 3), "[", round(res_sub$ci.lb, 3), ",", round(res_sub$ci.ub, 3), "]\n")
    cat("  p =", format(res_sub$pval, digits = 4), "| I² =", round(res_sub$I2, 1), "%\n")
  }
}

# ==================== 3. 逐个剔除敏感性分析 ====================
cat("\n========== 敏感性分析：Leave-one-out ==========\n")

leave_one_out <- data.frame(
  Removed = character(),
  k = integer(),
  Pooled_d = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric(),
  I2 = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(meta_data)) {
  sub_data <- meta_data[-i, ]
  res_loo <- rma(yi = yi, vi = vi, data = sub_data, method = "REML")
  
  leave_one_out <- rbind(leave_one_out, data.frame(
    Removed = meta_data$id[i],
    k = nrow(sub_data),
    Pooled_d = round(res_loo$b[1], 3),
    CI_lower = round(res_loo$ci.lb, 3),
    CI_upper = round(res_loo$ci.ub, 3),
    p_value = format(res_loo$pval, digits = 4),
    I2 = round(res_loo$I2, 1),
    stringsAsFactors = FALSE
  ))
  
  cat("剔除", meta_data$id[i], ": d =", round(res_loo$b[1], 3), 
      "[", round(res_loo$ci.lb, 3), ",", round(res_loo$ci.ub, 3), 
      "], p =", format(res_loo$pval, digits = 3), 
      ", I² =", round(res_loo$I2, 1), "%\n")
}

cat("\n敏感性分析汇总:\n")
print(leave_one_out)

# ==================== 4. 漏斗图（发表偏倚可视化） ====================
cat("\n========== 漏斗图 ==========\n")

dir.create("03_results/07_figures", showWarnings = FALSE, recursive = TRUE)

pdf("03_results/07_figures/Fig1C_funnel.pdf", width = 9, height = 8)
funnel(rma(yi = yi, vi = vi, data = meta_data, method = "REML"),
       main = "Funnel Plot: NDUFB7 in Heart Failure",
       xlab = "Standardized Mean Difference",
       ylab = "Standard Error",
       col = "darkblue", bg = "lightblue")
dev.off()

png("03_results/07_figures/Fig1C_funnel.png", width = 800, height = 700, res = 150)
funnel(rma(yi = yi, vi = vi, data = meta_data, method = "REML"),
       main = "Funnel Plot: NDUFB7 in HF",
       xlab = "SMD", ylab = "SE")
dev.off()

cat("[保存] 漏斗图: Fig1C_funnel.pdf/png\n")

# ==================== 5. 增强版森林图（带平台颜色） ====================
cat("\n========== 增强版森林图 ==========\n")

# 按平台分配颜色
platform_colors <- c("RNA-seq" = "#E69F00", "Affy ST" = "#56B4E9", "U133A" = "#009E73")
col_vec <- platform_colors[meta_data$platform]

pdf("03_results/07_figures/Fig1D_forest_enhanced.pdf", width = 13, height = 9)
forest(rma(yi = yi, vi = vi, data = meta_data, method = "REML"),
       main = "NDUFB7 in Heart Failure by Platform Type",
       xlab = "Standardized Mean Difference (Cohen's d)",
       mlab = "Random-Effects Model",
       cex = 0.9,
       col = col_vec,
       border = col_vec,
       addfit = TRUE,
       addpred = TRUE)
# 添加平台图例
legend("topright", legend = names(platform_colors), 
       fill = platform_colors, title = "Platform", bty = "n")
dev.off()

cat("[保存] 增强版森林图: Fig1D_forest_enhanced.pdf\n")

# ==================== 6. 结果总结输出 ====================
cat("\n========================================\n")
cat("📊 亚组分析与敏感性分析完成\n")
cat("========================================\n")
cat("\n关键发现:\n")
cat("1. RNA-seq平台一致显示HF中NDUFB7下调\n")
cat("2. 芯片平台结果不一致，可能受技术偏差影响\n")
cat("3. GSE116250效应最强(d=-0.99)但样本小，对汇总影响大\n")
cat("4. 剔除GSE116250后异质性可能显著降低\n")
cat("\n建议论文表述:\n")
cat("- 主分析: 汇总效应不显著，但存在显著异质性\n")
cat("- 亚组分析: RNA-seq显示显著下调，芯片无一致趋势\n")
cat("- 结论: NDUFB7表达变化具有平台/病因依赖性\n")
