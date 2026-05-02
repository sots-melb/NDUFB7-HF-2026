#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V145_GSE57338_Etiology")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V145_FIX2: GSE57338病因标签+NDUFB7+铁死亡")
message("========================================")

# === 1. 加载本地gene-level RDS（已确认存在，46M）===
rds_file <- "03_results/V133_GSE57338/GSE57338_gene_level.rds"
if (!file.exists(rds_file)) {
  # 广泛搜索
  rds_file <- list.files(c("01_data","03_results"), pattern = "GSE57338.*gene_level.*\\.rds$", 
                        full.names = TRUE, recursive = TRUE)[1]
}
if (is.null(rds_file) || !file.exists(rds_file)) stop("GSE57338 gene_level RDS not found")

obj <- readRDS(rds_file)
if (is.list(obj) && "exprs" %in% names(obj)) {
  exprs <- obj$exprs
} else if (is.matrix(obj) || is.data.frame(obj)) {
  exprs <- as.matrix(obj)
} else {
  stop("RDS structure unrecognized")
}
message("[PASS] Exprs loaded: ", nrow(exprs), " genes × ", ncol(exprs), " samples")

# === 2. 从series matrix解析表型（手动解析313样本）===
sm_file <- "01_data/01_raw_geo/GSE57338/GSE57338_series_matrix.txt.gz"
pheno <- data.frame(Sample = colnames(exprs), Disease = NA_character_, stringsAsFactors = FALSE)

if (file.exists(sm_file)) {
  lines <- readLines(sm_file, n = 300, warn = FALSE)
  # 提取!Sample_characteristics_ch1 和 !Sample_title
  char_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
  title_lines <- lines[grep("^!Sample_title", lines)]
  
  # 解析为向量
  parse_line <- function(l) {
    parts <- strsplit(l, "\t")[[1]]
    if (length(parts) > 1) return(parts[-1]) else return(NULL)
  }
  
  char_vals <- parse_line(char_lines[1])
  title_vals <- parse_line(title_lines[1])
  
  if (length(char_vals) == ncol(exprs)) {
    v <- tolower(char_vals)
    disease <- rep(NA_character_, length(v))
    disease[grep("ischem|isch|icm", v)] <- "ICM"
    disease[grep("dilated|dcm|idiopathic", v)] <- "DCM"
    disease[grep("non-failing|nf|control|healthy|normal", v)] <- "NF"
    pheno$Disease <- disease
    message("[PASS] Disease parsed from characteristics_ch1")
  } else if (length(title_vals) == ncol(exprs)) {
    v <- tolower(title_vals)
    disease <- rep(NA_character_, length(v))
    disease[grep("isch|icm", v)] <- "ICM"
    disease[grep("dilat|dcm", v)] <- "DCM"
    disease[grep("control|nf|normal", v)] <- "NF"
    pheno$Disease <- disease
    message("[PASS] Disease parsed from title")
  }
}

# 如果仍全NA，尝试从GSE57338_series_matrix.txt.gz的SOFT格式解析更多字段
if (all(is.na(pheno$Disease))) {
  # 尝试source_name_ch1
  src_lines <- lines[grep("^!Sample_source_name_ch1", lines)]
  src_vals <- parse_line(src_lines[1])
  if (length(src_vals) == ncol(exprs)) {
    v <- tolower(src_vals)
    disease <- rep(NA_character_, length(v))
    disease[grep("isch|icm", v)] <- "ICM"
    disease[grep("dilat|dcm", v)] <- "DCM"
    disease[grep("control|nf|normal", v)] <- "NF"
    pheno$Disease <- disease
    message("[PASS] Disease parsed from source_name_ch1")
  }
}

message("\n--- Disease distribution ---")
print(table(pheno$Disease, useNA = "ifany"))

# === 3. 提取NDUFB7 ===
if (!("NDUFB7" %in% rownames(exprs))) stop("NDUFB7 not in rownames")
ndufb7_expr <- as.numeric(exprs["NDUFB7", ])
names(ndufb7_expr) <- colnames(exprs)

# === 4. 加载铁死亡评分（从V141）===
ferro_file <- "03_results/V141_GSE57338_Bulk_Fix/V141_bulk_ferroptosis.csv"
ferro_score <- NULL
if (file.exists(ferro_file)) {
  ferro_df <- read.csv(ferro_file, stringsAsFactors = FALSE)
  if (all(c("sample","ferroptosis_score") %in% colnames(ferro_df))) {
    ferro_score <- ferro_df$ferroptosis_score
    names(ferro_score) <- ferro_df$sample
    message("[PASS] Ferroptosis score loaded from V141")
  }
}

# === 5. 对齐并统计 ===
common <- intersect(names(ndufb7_expr), pheno$Sample)
df <- data.frame(
  Sample = common,
  Disease = factor(pheno$Disease[match(common, pheno$Sample)], levels = c("NF", "DCM", "ICM")),
  NDUFB7 = ndufb7_expr[common],
  stringsAsFactors = FALSE
)

if (!is.null(ferro_score)) {
  common2 <- intersect(common, names(ferro_score))
  df <- df[df$Sample %in% common2, ]
  df$Ferroptosis_Score <- ferro_score[df$Sample]
}

valid <- !is.na(df$Disease) & is.finite(df$NDUFB7)
df <- df[valid, ]

message("\n--- Final analysis set ---")
print(table(df$Disease))

if (length(unique(df$Disease)) >= 2 && nrow(df) > 10) {
  kw <- kruskal.test(NDUFB7 ~ Disease, data = df)
  message("\nKruskal-Wallis: χ² = ", round(kw$statistic, 2), 
          ", df = ", kw$parameter,
          ", p = ", format(kw$p.value, digits=2, scientific=TRUE))
  
  if (length(unique(df$Disease)) >= 3) {
    pw <- pairwise.wilcox.test(df$NDUFB7, df$Disease, p.adjust.method = "BH", exact = FALSE)
    message("\nPairwise Wilcoxon (BH):")
    print(pw$p.value)
    write.csv(as.data.frame(pw$p.value), file.path(outdir, "V145_pairwise_pvalues.csv"))
  }
  
  # 效应量
  if (all(c("NF","ICM") %in% df$Disease)) {
    nf <- df$NDUFB7[df$Disease == "NF"]
    icm <- df$NDUFB7[df$Disease == "ICM"]
    cohen_d <- (mean(nf, na.rm=TRUE) - mean(icm, na.rm=TRUE)) / 
               sqrt(((length(nf)-1)*var(nf, na.rm=TRUE) + (length(icm)-1)*var(icm, na.rm=TRUE)) / (length(nf)+length(icm)-2))
    message("\nNF vs ICM Cohen's d: ", round(cohen_d, 3))
  }
  
  write.csv(df, file.path(outdir, "V145_GSE57338_NDUFB7_by_Etiology.csv"), row.names = FALSE)
  
  # 投稿级Boxplot
  p <- ggplot(df, aes(x = Disease, y = NDUFB7, fill = Disease)) +
    geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.12, alpha = 0.35, size = 1.2, color = "black") +
    scale_fill_manual(values = c("NF" = "#0072B2", "DCM" = "#E69F00", "ICM" = "#D55E00"),
                      labels = c("NF" = "Non-Failing", "DCM" = "Dilated CM", "ICM" = "Ischemic CM")) +
    labs(
      title = "NDUFB7 Expression by Etiology",
      subtitle = paste0("GSE57338 | Kruskal-Wallis p = ", signif(kw$p.value, 2), " | n = ", nrow(df)),
      x = "Etiology", y = "NDUFB7 Expression (log2)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      panel.background = element_rect(fill = "#F5F5F5", color = NA),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(face = "bold")
    )
  
  ggsave(file.path(outdir, "V145_Etiology_Boxplot_300dpi.png"), p, width = 6, height = 5, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "V145_Etiology_Boxplot.pdf"), p, width = 6, height = 5, device = cairo_pdf)
  
  # 铁死亡相关性（如果可用）
  if ("Ferroptosis_Score" %in% colnames(df)) {
    ct <- cor.test(df$NDUFB7, df$Ferroptosis_Score, method = "spearman", use = "complete.obs")
    message("\nNDUFB7 vs Ferroptosis (by etiology): rho=", round(ct$estimate, 3), 
            " p=", format(ct$p.value, digits=2, scientific=TRUE))
  }
  
  message("[DONE] V145_FIX2: ", outdir)
} else {
  message("[FAIL] Insufficient groups or samples for comparison")
  message("  Unique diseases: ", paste(unique(df$Disease), collapse = ", "))
  message("  Total samples: ", nrow(df))
}
