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
message("V145_FIX: GSE57338病因标签精确分组")
message("========================================")

# ========== 1. 加载GSE57338（本地优先）==========
gse <- NULL
sm_local <- file.path(PROJECT, "01_data/01_raw_geo/GSE57338/GSE57338_series_matrix.txt.gz")

if (file.exists(sm_local)) {
  message("[LOCAL] Loading from: ", sm_local)
  tryCatch({
    gse_list <- getGEO(filename = sm_local, getGPL = FALSE)
    if (is.list(gse_list) && length(gse_list) > 0) {
      gse <- gse_list[[1]]
      message("[PASS] Local series matrix loaded")
    }
  }, error = function(e) message("[WARN] Local parse failed: ", conditionMessage(e)))
}

if (is.null(gse)) {
  message("[ONLINE] Attempting GEO online...")
  tryCatch({
    gse_list <- getGEO("GSE57338", GSEMatrix = TRUE, getGPL = FALSE)
    if (is.list(gse_list) && length(gse_list) > 0 && inherits(gse_list[[1]], "ExpressionSet")) {
      gse <- gse_list[[1]]
      message("[PASS] Online GEO loaded")
    } else {
      message("[WARN] Online returned non-ExpressionSet: ", class(gse_list))
    }
  }, error = function(e) message("[FAIL] Online: ", conditionMessage(e)))
}

if (is.null(gse) || !inherits(gse, "ExpressionSet")) {
  # 最终回退：手动解析series matrix表头
  message("[FALLBACK] Manual parsing series matrix header...")
  if (file.exists(sm_local)) {
    lines <- readLines(sm_local, n = 200, warn = FALSE)
    # 提取样本信息行
    sample_lines <- lines[grep("^!Sample_", lines)]
    # 解析为数据框
    parse_geo_header <- function(lines) {
      vals <- strsplit(lines, "\t")
      keys <- sapply(vals, function(x) x[1])
      # 去重键，保留第一个
      unique_keys <- unique(keys)
      df <- data.frame(row.names = 1:(length(vals[[1]])-1), stringsAsFactors = FALSE)
      for (k in unique_keys) {
        idx <- which(keys == k)[1]
        row_vals <- vals[[idx]][-1]
        col_name <- gsub("^!Sample_", "", k)
        if (length(row_vals) > 0) {
          df[[col_name]] <- row_vals
        }
      }
      return(df)
    }
    pheno <- parse_geo_header(sample_lines)
    message("[PASS] Manual parse: ", nrow(pheno), " samples")
    
    # 加载exprs（使用已有的RDS）
    rds_file <- "03_results/V133_GSE57338/GSE57338_gene_level.rds"
    if (file.exists(rds_file)) {
      obj <- readRDS(rds_file)
      if (is.list(obj) && "exprs" %in% names(obj)) {
        exprs_mat <- obj$exprs
      } else {
        exprs_mat <- as.matrix(obj)
      }
      message("[PASS] Exprs from RDS: ", ncol(exprs_mat), " samples")
    } else {
      stop("No expression matrix available")
    }
  } else {
    stop("GSE57338 series matrix not found locally and online failed")
  }
} else {
  pheno <- pData(gse)
  exprs_mat <- exprs(gse)
  message("[PASS] Standard GEO: ", nrow(pheno), " pheno × ", ncol(exprs_mat), " exprs")
}

# ========== 2. 解析病因标签（多字段尝试）==========
disease_label <- rep(NA_character_, nrow(pheno))

# 尝试1: characteristics列
char_cols <- grep("characteristics", colnames(pheno), value = TRUE, ignore.case = TRUE)
for (cc in char_cols) {
  v <- tolower(as.character(pheno[[cc]]))
  disease_label[grep("ischem|isch|icm", v)] <- "ICM"
  disease_label[grep("dilated|dcm|dilated cardiomyopathy", v)] <- "DCM"
  disease_label[grep("non-failing|non failing|nf|control|healthy|normal", v)] <- "NF"
}

# 尝试2: title列
if (all(is.na(disease_label))) {
  t <- tolower(as.character(pheno$title))
  disease_label[grep("isch|icm", t)] <- "ICM"
  disease_label[grep("dilat|dcm", t)] <- "DCM"
  disease_label[grep("control|nf|normal", t)] <- "NF"
}

# 尝试3: source_name_ch1
if (all(is.na(disease_label))) {
  s <- tolower(as.character(pheno$source_name_ch1))
  disease_label[grep("isch|icm", s)] <- "ICM"
  disease_label[grep("dilat|dcm", s)] <- "DCM"
  disease_label[grep("control|nf|normal", s)] <- "NF"
}

# 尝试4: 任何包含disease的列
if (all(is.na(disease_label))) {
  dis_cols <- grep("disease|diagnosis|etiology|group", colnames(pheno), value = TRUE, ignore.case = TRUE)
  for (cc in dis_cols) {
    v <- tolower(as.character(pheno[[cc]]))
    disease_label[grep("ischem|isch|icm", v)] <- "ICM"
    disease_label[grep("dilated|dcm", v)] <- "DCM"
    disease_label[grep("non-fail|nf|control|healthy", v)] <- "NF"
  }
}

message("\n--- Disease label distribution ---")
print(table(disease_label, useNA = "ifany"))

# 调试：如果全NA，输出前5行表型供诊断
if (all(is.na(disease_label))) {
  message("\n[DEBUG] First 5 rows of pheno:")
  print(head(pheno[, 1:min(5, ncol(pheno))]))
  message("\n[DEBUG] All column names:")
  print(colnames(pheno))
}

# ========== 3. 提取NDUFB7 ==========
ndufb7_expr <- NULL
if ("NDUFB7" %in% rownames(exprs_mat)) {
  ndufb7_expr <- as.numeric(exprs_mat["NDUFB7", ])
  names(ndufb7_expr) <- colnames(exprs_mat)
  message("[PASS] NDUFB7 from rownames")
} else {
  # 探针映射回退
  pmap_file <- "03_results/V140_BIC_Resolution/V133_NDUFB7_probes.csv"
  if (file.exists(pmap_file)) {
    pmap <- fread(pmap_file)
    probes <- pmap$ID
    rows <- rownames(exprs_mat) %in% probes
    if (any(rows)) {
      ndufb7_expr <- colMeans(exprs_mat[rows, , drop = FALSE], na.rm = TRUE)
      names(ndufb7_expr) <- colnames(exprs_mat)
      message("[INFO] NDUFB7 via probe map")
    }
  }
}

if (is.null(ndufb7_expr)) {
  hits <- grep("NDUFB7", rownames(exprs_mat), ignore.case = TRUE)
  if (length(hits) > 0) {
    ndufb7_expr <- as.numeric(exprs_mat[hits[1], ])
    names(ndufb7_expr) <- colnames(exprs_mat)
    message("[WARN] Fuzzy match: ", rownames(exprs_mat)[hits[1]])
  }
}

if (is.null(ndufb7_expr)) stop("NDUFB7 extraction failed")

# ========== 4. 对齐样本 ==========
# pheno行名通常是GSMxxx
pheno_names <- rownames(pheno)
exprs_names <- names(ndufb7_expr)

# 如果exprs_names是NULL，用位置对齐
if (is.null(exprs_names) || length(exprs_names) == 0) {
  if (nrow(pheno) == length(ndufb7_expr)) {
    exprs_names <- pheno_names
    names(ndufb7_expr) <- exprs_names
    message("[WARN] Positional alignment (no names in exprs)")
  }
}

common <- intersect(pheno_names, exprs_names)
if (length(common) == 0) {
  # 尝试通过title匹配
  if ("title" %in% colnames(pheno)) {
    common <- intersect(pheno$title, exprs_names)
  }
}
if (length(common) == 0 && nrow(pheno) == length(ndufb7_expr)) {
  common <- pheno_names
  names(ndufb7_expr) <- common
  message("[WARN] Forced positional alignment")
}

disease_aligned <- disease_label[match(common, pheno_names)]
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

# ========== 5. 统计检验与可视化 ==========
if (length(unique(df$Disease)) >= 2 && nrow(df) > 10) {
  kw <- kruskal.test(NDUFB7 ~ Disease, data = df)
  message("\nKruskal-Wallis: χ² = ", round(kw$statistic, 2), 
          ", df = ", kw$parameter,
          ", p = ", format(kw$p.value, digits=2, scientific=TRUE))
  
  # 两两Wilcoxon + BH
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
  
  fwrite(df, file.path(outdir, "V145_GSE57338_NDUFB7_by_Etiology.csv"))
  
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
  
  message("[DONE] V145_FIX: ", outdir)
} else {
  message("[FAIL] Insufficient groups or samples")
  message("  Unique diseases: ", paste(unique(df$Disease), collapse = ", "))
  message("  Total samples: ", nrow(df))
}
