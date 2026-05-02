#!/usr/bin/env Rscript
# V86_T26: NDUFAF3-NDUFB7共表达验证
# 支撑"组装缺陷轴"假说

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

# --- 查找GSE183852数据 ---
files <- list.files("Downloads", pattern = "GSE183852.*\\.(rds|Robj)", full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) {
  files <- list.files("01_data", pattern = "GSE183852", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
}
if (length(files) == 0) stop("[ERROR] 未找到GSE183852数据文件")

srt <- readRDS(files[1])
message("=== T26: NDUFAF3-NDUFB7共表达验证 ===")
message("加载: ", basename(files[1]), " | ", ncol(srt), " cells/nuclei")

# --- 提取CM子集 ---
if ("cell_type" %in% colnames(srt@meta.data)) {
  cm <- subset(srt, cell_type %in% c("Cardiomyocyte", "CM"))
} else if ("predicted.id" %in% colnames(srt@meta.data)) {
  cm <- subset(srt, predicted.id == "Cardiomyocyte")
} else {
  cm <- srt
  message("[WARN] 无细胞类型注释，使用全部细胞")
}
message("CM子集: ", ncol(cm), " cells")

# --- 提取基因表达 ---
genes <- c("NDUFAF3", "NDUFB7")
avail <- intersect(genes, rownames(cm))

if (length(avail) < 2) {
  message("[FAIL] 缺失基因。可用: ", paste(avail, collapse = ", "))
  all_g <- rownames(cm)
  for (g in genes) {
    m <- grep(paste0("^", g, "$"), all_g, ignore.case = TRUE, value = TRUE)
    if (length(m) == 0) m <- grep(g, all_g, ignore.case = TRUE, value = TRUE)
    message("  可能的匹配 '", g, "': ", paste(m, collapse = ", "))
  }
  stop("基因名不匹配，请检查")
}

expr <- FetchData(cm, vars = avail)
cor_pearson <- cor.test(expr[[1]], expr[[2]], method = "pearson")
cor_spear <- cor.test(expr[[1]], expr[[2]], method = "spearman")

message("")
message("=== T26 统计结果 ===")
message("Pearson r  = ", round(cor_pearson$estimate, 3), "  p = ", format(cor_pearson$p.value, digits = 2, scientific = TRUE))
message("Spearman ρ = ", round(cor_spear$estimate, 3),     "  p = ", format(cor_spear$p.value, digits = 2, scientific = TRUE))
message("N cells  = ", nrow(expr))

# --- 可视化 ---
p <- ggplot(expr, aes_string(x = avail[1], y = avail[2])) +
  geom_point(alpha = 0.3, size = 0.5, color = "#31688E") +
  geom_smooth(method = "lm", color = "#FDE725", se = TRUE, linewidth = 1) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("r = ", round(cor_pearson$estimate, 3), 
                         "\np = ", format(cor_pearson$p.value, digits = 1, scientific = TRUE)),
           hjust = 1.1, vjust = -0.5, size = 3) +
  labs(title = "NDUFAF3-NDUFB7 Co-expression in Cardiomyocytes",
       x = paste0(avail[1], " expression (log-norm)"),
       y = paste0(avail[2], " expression (log-norm)")) +
  theme_minimal(base_size = 10)

outdir <- file.path(PROJECT_DIR, "03_results/T26_NDUFAF3_NDUFB7")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(outdir, "T26_scatter.png"), p, width = 6, height = 5, dpi = 300)

stats <- data.frame(
  Metric = c("Pearson_r", "Pearson_p", "Spearman_rho", "Spearman_p", "N_cells"),
  Value = c(cor_pearson$estimate, cor_pearson$p.value, 
            cor_spear$estimate, cor_spear$p.value, nrow(expr))
)
write.csv(stats, file.path(outdir, "T26_stats.csv"), row.names = FALSE)

# --- 判断 ---
message("")
if (cor_pearson$estimate > 0.6 && cor_pearson$p.value < 0.05) {
  message("[PASS] 共表达显著且强相关！直接支持'组装缺陷轴'假说")
  message("[IMPACT] 可写入Results，支撑Fig 4B")
} else if (cor_pearson$p.value < 0.05) {
  message("[PARTIAL] 共表达存在但强度中等，可写入但需保守表述")
} else {
  message("[WARN] 共表达不显著，建议启动替代方案:")
  message("  方案A: 检查NDUFAF3基因名/Ensembl ID匹配")
  message("  方案B: 换用GSE57338 Bulk数据验证")
  message("  方案C: 从叙事中弱化NDUFAF3轴，聚焦铁死亡")
}
message("[DONE] 结果保存: ", outdir)
