#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(ggplot2)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V145_GSE57338_Etiology")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V145: GSE57338病因标签精确分组")
message("========================================")

# ========== 1. 加载GSE57338表型 ==========
gse <- tryCatch({
  getGEO("GSE57338", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
}, error = function(e) {
  sm <- file.path(PROJECT, "01_data/01_raw_geo/GSE57338/GSE57338_series_matrix.txt.gz")
  if (file.exists(sm)) getGEO(filename = sm, getGPL = FALSE)[[1]] else NULL
})

if (is.null(gse)) stop("GSE57338 unavailable")

pheno <- pData(gse)
exprs_mat <- exprs(gse)

message("[PASS] Pheno: ", nrow(pheno), " samples | Exprs: ", ncol(exprs_mat), " samples")

# ========== 2. 解析病因标签 ==========
# GSE57338常见标签在characteristics_ch1/ch2中
disease_label <- rep(NA_character_, nrow(pheno))

# 搜索所有characteristics列
char_cols <- grep("characteristics", colnames(pheno), value = TRUE)
for (cc in char_cols) {
  v <- tolower(as.character(pheno[[cc]]))
  # ICM模式
  disease_label[grep("ischemic|icm|isch", v)] <- "ICM"
  # DCM模式
  disease_label[grep("dilated|dcm|dilated cardiomyopathy", v)] <- "DCM"
  # NF/Control模式
  disease_label[grep("non-failing|nf|control|healthy|normal", v)] <- "NF"
}

# 如果characteristics未解析成功，尝试title列
if (all(is.na(disease_label))) {
  t <- tolower(as.character(pheno$title))
  disease_label[grep("isch|icm", t)] <- "ICM"
  disease_label[grep("dilat|dcm", t)] <- "DCM"
  disease_label[grep("control|nf|normal", t)] <- "NF"
}

# 如果仍然失败，尝试source_name_ch1
if (all(is.na(disease_label))) {
  s <- tolower(as.character(pheno$source_name_ch1))
  disease_label[grep("isch|icm", s)] <- "ICM"
  disease_label[grep("dilat|dcm", s)] <- "DCM"
  disease_label[grep("control|nf|normal", s)] <- "NF"
}

# 样本级调试输出
message("\n--- Disease label distribution ---")
print(table(disease_label, useNA = "ifany"))

# ========== 3. 提取NDUFB7 ==========
# 需要探针映射（GSE57338使用GPL11532）
# 简化：如果V133有probe map，使用；否则用第一个含NDUFB7的探针
ndufb7_expr <- NULL
if ("NDUFB7" %in% rownames(exprs_mat)) {
  ndufb7_expr <- as.numeric(exprs_mat["NDUFB7", ])
} else {
  # 尝试从V133加载映射
  pmap_file <- "03_results/V140_BIC_Resolution/V133_NDUFB7_probes.csv"
  if (file.exists(pmap_file)) {
    pmap <- fread(pmap_file)
    probes <- pmap$ID
    rows <- rownames(exprs_mat) %in% probes
    if (any(rows)) {
      ndufb7_expr <- colMeans(exprs_mat[rows, , drop = FALSE], na.rm = TRUE)
      message("[INFO] NDUFB7 extracted via probe map")
    }
  }
}

if (is.null(ndufb7_expr)) {
  # 最终回退：行名模糊匹配
  hits <- grep("NDUFB7", rownames(exprs_mat), ignore.case = TRUE)
  if (length(hits) > 0) {
    ndufb7_expr <- as.numeric(exprs_mat[hits[1], ])
    message("[WARN] NDUFB7 matched by fuzzy search: ", rownames(exprs_mat)[hits[1]])
  }
}

if (is.null(ndufb7_expr)) stop("NDUFB7 not found in expression matrix")

# ========== 4. 对齐样本 ==========
# pheno行名通常是GSMxxx，exprs列名也应该是GSMxxx
common <- intersect(rownames(pheno), colnames(exprs_mat))
if (length(common) == 0) {
  # 尝试通过title或其他字段匹配
  common <- intersect(rownames(pheno), names(ndufb7_expr))
}
if (length(common) == 0) {
  # 位置对齐回退
  if (nrow(pheno) == length(ndufb7_expr)) {
    common <- rownames(pheno)
    names(ndufb7_expr) <- common
    message("[WARN] Positional alignment used")
  }
}

disease_aligned <- disease_label[match(common, rownames(pheno))]
ndufb7_aligned <- ndufb7_expr[common]

valid <- !is.na(disease_aligned) & is.finite(ndufb7_aligned)
df <- data.frame(
  Sample = common[valid],
  Disease = factor(disease_aligned[valid], levels = c("NF", "DCM", "ICM")),
  NDUFB7 = ndufb7_aligned[valid],
  stringsAsFactors = FALSE
)

message("\n--- Final analysis set ---")
print(table(df$Disease))

# ========== 5. 统计检验 ==========
if (length(unique(df$Disease)) >= 2) {
  # Kruskal-Wallis
  kw <- kruskal.test(NDUFB7 ~ Disease, data = df)
  message("\nKruskal-Wallis: χ² = ", round(kw$statistic, 2), 
          ", df = ", kw$parameter,
          ", p = ", format(kw$p.value, digits=2, scientific=TRUE))
  
  # 两两比较（Wilcoxon + BH校正）
  pw <- pairwise.wilcox.test(df$NDUFB7, df$Disease, p.adjust.method = "BH", exact = FALSE)
  message("\nPairwise Wilcoxon (BH-adjusted):")
  print(pw$p.value)
  
  # 效应量：NF vs ICM
  if (all(c("NF","ICM") %in% df$Disease)) {
    nf <- df$NDUFB7[df$Disease == "NF"]
    icm <- df$NDUFB7[df$Disease == "ICM"]
    cohen_d <- (mean(nf, na.rm=TRUE) - mean(icm, na.rm=TRUE)) / 
               sqrt(((length(nf)-1)*var(nf, na.rm=TRUE) + (length(icm)-1)*var(icm, na.rm=TRUE)) / (length(nf)+length(icm)-2))
    message("\nNF vs ICM Cohen's d: ", round(cohen_d, 3))
  }
  
  # 保存
  fwrite(df, file.path(outdir, "V145_GSE57338_NDUFB7_by_Etiology.csv"))
  write.csv(as.data.frame(pw$p.value), file.path(outdir, "V145_pairwise_pvalues.csv"))
  
  # 三分类Boxplot
  p <- ggplot(df, aes(x = Disease, y = NDUFB7, fill = Disease)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
    scale_fill_manual(values = c("NF" = "#0072B2", "DCM" = "#E69F00", "ICM" = "#D55E00")) +
    labs(
      title = "NDUFB7 Expression by Etiology (GSE57338)",
      subtitle = paste0("Kruskal-Wallis p = ", signif(kw$p.value, 2)),
      x = "Etiology", y = "NDUFB7 Expression"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      panel.background = element_rect(fill = "#F0F0F0", color = NA),
      panel.grid.major = element_line(color = "white")
    )
  
  ggsave(file.path(outdir, "V145_Etiology_Boxplot.png"), p, width = 6, height = 5, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "V145_Etiology_Boxplot.pdf"), p, width = 6, height = 5, device = cairo_pdf)
  
  message("[DONE] V145: ", outdir)
} else {
  message("[FAIL] Insufficient disease groups for comparison")
}
