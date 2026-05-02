#!/usr/bin/env Rscript
# V88: 诊断GSE183852结构 + 自适应加载 + T26共表达验证

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Biobase)  # for ExpressionSet
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

message("========================================")
message("V88: GSE183852 结构诊断 + T26修复")
message("========================================")

# --- 1. 查找所有GSE183852相关文件 ---
files <- c(
  list.files("Downloads", pattern = "GSE183852", full.names = TRUE, ignore.case = TRUE, recursive = TRUE),
  list.files("01_data", pattern = "GSE183852", full.names = TRUE, ignore.case = TRUE, recursive = TRUE),
  list.files("03_results", pattern = "GSE183852", full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
)
files <- unique(files[file.info(files)$size > 1000])  # 排除空文件
message("找到 ", length(files), " 个GSE183852候选文件:")
for (f in files) message("  ", basename(f), " (", round(file.info(f)$size/1024/1024, 1), " MB)")

if (length(files) == 0) stop("[FAIL] 未找到任何GSE183852文件")

# --- 2. 逐个探测结构 ---
for (f in files) {
  message("")
  message(">>> 探测: ", basename(f))
  obj <- readRDS(f)
  cls <- class(obj)
  message("  对象类型: ", paste(cls, collapse = ", "))
  
  if ("Seurat" %in% cls) {
    message("  [FOUND] 标准Seurat对象，细胞数: ", ncol(obj))
    BEST_FILE <- f
    BEST_TYPE <- "Seurat"
    break
  }
  
  if ("ExpressionSet" %in% cls) {
    message("  [FOUND] ExpressionSet (Bulk), 样本数: ", ncol(obj), ", 基因数: ", nrow(obj))
    BEST_FILE <- f
    BEST_TYPE <- "ExpressionSet"
    break
  }
  
  if (is.list(obj)) {
    message("  [LIST] List结构，长度: ", length(obj))
    message("  名称: ", paste(names(obj)[1:min(5, length(obj))], collapse = ", "))
    
    # 检查是否包含Seurat
    if (any(sapply(obj, function(x) "Seurat" %in% class(x)))) {
      idx <- which(sapply(obj, function(x) "Seurat" %in% class(x)))[1]
      message("  [FOUND] List中包含Seurat对象，取第", idx, "个")
      obj <- obj[[idx]]
      BEST_FILE <- f
      BEST_TYPE <- "Seurat_in_list"
      break
    }
    
    # 检查是否包含ExpressionSet
    if (any(sapply(obj, function(x) "ExpressionSet" %in% class(x)))) {
      idx <- which(sapply(obj, function(x) "ExpressionSet" %in% class(x)))[1]
      message("  [FOUND] List中包含ExpressionSet，取第", idx, "个")
      obj <- obj[[idx]]
      BEST_FILE <- f
      BEST_TYPE <- "ExpressionSet_in_list"
      break
    }
    
    # 检查是否有counts和meta
    if (all(c("counts", "metadata") %in% names(obj))) {
      message("  [FOUND] List包含counts+metadata，可手动构建Seurat")
      BEST_FILE <- f
      BEST_TYPE <- "counts_meta_list"
      break
    }
  }
}

if (!exists("BEST_TYPE")) stop("[FAIL] 未识别可用数据结构")

message("")
message("========================================")
message("最终选择: ", basename(BEST_FILE))
message("数据类型: ", BEST_TYPE)
message("========================================")

# --- 3. 统一转换为Seurat对象 ---
if (BEST_TYPE == "Seurat" || BEST_TYPE == "Seurat_in_list") {
  srt <- obj
  message("[PASS] 直接使用Seurat对象")
  
} else if (BEST_TYPE %in% c("ExpressionSet", "ExpressionSet_in_list")) {
  # Bulk ExpressionSet -> 创建伪单细胞Seurat（用于基因表达提取）
  exprs_mat <- exprs(obj)  # 基因×样本矩阵
  pheno <- pData(obj)
  
  message("[CONVERT] ExpressionSet -> Seurat (", nrow(exprs_mat), " genes × ", ncol(exprs_mat), " samples)")
  
  srt <- CreateSeuratObject(counts = exprs_mat, meta.data = pheno)
  # 标准化
  srt <- NormalizeData(srt)
  
} else if (BEST_TYPE == "counts_meta_list") {
  counts <- obj$counts
  meta <- obj$metadata
  srt <- CreateSeuratObject(counts = counts, meta.data = meta)
  srt <- NormalizeData(srt)
  message("[CONVERT] List(counts+meta) -> Seurat")
}

message("Seurat对象: ", nrow(srt), " features × ", ncol(srt), " cells/samples")

# --- 4. 提取CM子集（自适应） ---
meta_cols <- colnames(srt@meta.data)
message("")
message("Meta.data列: ", paste(meta_cols, collapse = ", "))

cm_keywords <- c("cell_type", "predicted.id", "celltype", "CellType", "cell_type1", "singleR.labels", "seurat_clusters")
cm_col <- intersect(cm_keywords, meta_cols)[1]

if (!is.na(cm_col)) {
  message("[INFO] 使用列 '", cm_col, "' 识别CM")
  cm_vals <- unique(srt@meta.data[[cm_col]])
  message("  唯一值: ", paste(cm_vals, collapse = ", "))
  
  cm_patterns <- c("Cardiomyocyte", "CM", "cardio", "Cardio", "myocyte")
  cm_mask <- sapply(srt@meta.data[[cm_col]], function(x) any(grepl(paste(cm_patterns, collapse = "|"), as.character(x), ignore.case = TRUE)))
  
  if (sum(cm_mask) > 0) {
    cm <- subset(srt, cells = colnames(srt)[cm_mask])
    message("[PASS] CM子集: ", ncol(cm), " cells/samples")
  } else {
    message("[WARN] 未匹配到CM关键词，使用全部数据")
    cm <- srt
  }
} else {
  message("[WARN] 无细胞类型列，使用全部数据")
  cm <- srt
}

# --- 5. T26: NDUFAF3-NDUFB7共表达分析 ---
message("")
message(">>> T26: NDUFAF3-NDUFB7共表达验证")

genes <- c("NDUFAF3", "NDUFB7")
avail <- intersect(genes, rownames(cm))

if (length(avail) < 2) {
  message("[FAIL] 缺失基因。可用: ", paste(avail, collapse = ", "))
  all_g <- rownames(cm)
  for (g in genes) {
    m <- grep(paste0("^", g, "$"), all_g, ignore.case = TRUE, value = TRUE)
    if (length(m) == 0) m <- grep(g, all_g, ignore.case = TRUE, value = TRUE)
    message("  可能的匹配 '", g, "': ", paste(head(m, 5), collapse = ", "))
  }
  stop("基因名不匹配")
}

# 提取表达（Seurat v4/v5兼容）
expr <- FetchData(cm, vars = avail)
colnames(expr) <- avail

cor_pearson <- cor.test(expr[[1]], expr[[2]], method = "pearson")
cor_spear <- cor.test(expr[[1]], expr[[2]], method = "spearman")

message("")
message("=== T26 统计结果 ===")
message("Pearson r  = ", round(cor_pearson$estimate, 3), "  p = ", format(cor_pearson$p.value, digits = 2, scientific = TRUE))
message("Spearman ρ = ", round(cor_spear$estimate, 3),     "  p = ", format(cor_spear$p.value, digits = 2, scientific = TRUE))
message("N = ", nrow(expr))

# --- 6. 可视化 ---
p <- ggplot(expr, aes(x = .data[[avail[1]]], y = .data[[avail[2]]])) +
  geom_point(alpha = 0.3, size = 0.5, color = "#31688E") +
  geom_smooth(method = "lm", color = "#FDE725", se = TRUE, linewidth = 1) +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("r = ", round(cor_pearson$estimate, 3), 
                         "\np = ", format(cor_pearson$p.value, digits = 1, scientific = TRUE)),
           hjust = 1.1, vjust = -0.5, size = 3) +
  labs(title = "NDUFAF3-NDUFB7 Co-expression",
       subtitle = paste0("N = ", nrow(expr), " | Data: ", basename(BEST_FILE)),
       x = paste0(avail[1], " expression"),
       y = paste0(avail[2], " expression")) +
  theme_minimal(base_size = 10)

outdir <- file.path(PROJECT_DIR, "03_results/T26_NDUFAF3_NDUFB7")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(outdir, "V88_T26_scatter.png"), p, width = 6, height = 5, dpi = 300)

stats <- data.frame(
  File = basename(BEST_FILE),
  Data_Type = BEST_TYPE,
  N = nrow(expr),
  Pearson_r = cor_pearson$estimate,
  Pearson_p = cor_pearson$p.value,
  Spearman_rho = cor_spear$estimate,
  Spearman_p = cor_spear$p.value
)
write.csv(stats, file.path(outdir, "V88_T26_stats.csv"), row.names = FALSE)

# --- 7. 判断 ---
message("")
if (cor_pearson$estimate > 0.6 && cor_pearson$p.value < 0.05) {
  message("[PASS] 共表达显著且强相关！直接支持'组装缺陷轴'假说")
} else if (cor_pearson$p.value < 0.05) {
  message("[PARTIAL] 共表达存在但强度中等，可写入但需保守表述")
} else {
  message("[WARN] 共表达不显著")
}
message("[DONE] 结果保存: ", outdir)
