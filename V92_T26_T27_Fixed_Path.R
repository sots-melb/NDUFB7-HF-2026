#!/usr/bin/env Rscript
# V92: 固定路径执行T26/T27，禁止任何外部搜索或自适应加载

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(reshape2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

# --- 严格固定路径 ---
STD_DIR <- "01_data/02_single_cell/GSE183852"
FILE_P0 <- file.path(STD_DIR, "GSE183852_CM_annotated.rds")      # 49M, V81结果, 首选
FILE_P1 <- file.path(STD_DIR, "GSE183852_Pure_CM_Subsampled.RDS") # 915M, 备选
FILE_P2 <- file.path(STD_DIR, "GSE183852_scored.RDS")             # 8.3G, 次选
FILE_P3 <- file.path(STD_DIR, "GSE183852_DCM_Nuclei.Robj.gz")     # 8.2G, 末选

message("========================================")
message("V92: 固定路径 T26/T27")
message("========================================")

# --- 1. 按优先级加载 ---
FILE <- NULL
LABEL <- NULL

if (file.exists(FILE_P0)) {
  FILE <- FILE_P0; LABEL <- "P0_CM_annotated (V81结果)"
} else if (file.exists(FILE_P1)) {
  FILE <- FILE_P1; LABEL <- "P1_Pure_CM_Subsampled"
} else if (file.exists(FILE_P2)) {
  FILE <- FILE_P2; LABEL <- "P2_scored"
} else if (file.exists(FILE_P3)) {
  FILE <- FILE_P3; LABEL <- "P3_DCM_Nuclei_raw"
}

if (is.null(FILE)) {
  message("[FAIL] 标准目录内无可用文件: ", STD_DIR)
  message("标准目录应包含以下至少一项:")
  message("  - GSE183852_CM_annotated.rds")
  message("  - GSE183852_Pure_CM_Subsampled.RDS")
  message("  - GSE183852_scored.RDS")
  message("  - GSE183852_DCM_Nuclei.Robj.gz")
  message("[ACTION] 请先执行 V91 归档脚本")
  stop("数据未归档")
}

message("[LOAD] 数据源: ", LABEL)
message("[LOAD] 路径: ", FILE)
message("[LOAD] 读取中...")

# --- 2. 读取（按格式）---
if (grepl("\\.gz$", FILE, ignore.case = TRUE)) {
  con <- gzfile(FILE, "rb")
  obj <- readRDS(con)
  close(con)
} else {
  obj <- readRDS(FILE)
}

cls <- class(obj)
if ("Seurat" %in% cls) {
  srt <- obj
  message("[PASS] Seurat对象，", ncol(srt), " cells × ", nrow(srt), " genes")
} else if (is.list(obj) && any(sapply(obj, function(x) "Seurat" %in% class(x)))) {
  idx <- which(sapply(obj, function(x) "Seurat" %in% class(x)))[1]
  srt <- obj[[idx]]
  message("[PASS] List包裹Seurat，提取[", idx, "]，", ncol(srt), " cells")
} else {
  stop("[FAIL] 非Seurat对象: ", paste(cls, collapse = ", "))
}

# --- 3. CM确认 ---
# P0/P1 文件名已提示是CM子集，跳过提取
if (grepl("CM_annotated|Pure_CM|Subsampled", LABEL, ignore.case = TRUE)) {
  cm <- srt
  message("[INFO] 数据源已是CM子集，直接使用: ", ncol(cm), " cells")
} else {
  # P2/P3 需提取CM
  meta_cols <- colnames(srt@meta.data)
  cm_col <- intersect(c("cell_type","predicted.id","celltype"), meta_cols)[1]
  if (!is.na(cm_col)) {
    mask <- grepl("Cardio|CM|myocyte", srt@meta.data[[cm_col]], ignore.case = TRUE)
    if (sum(mask) > 0) {
      cm <- subset(srt, cells = colnames(srt)[mask])
      message("[SUBSET] CM提取: ", ncol(cm), "/", ncol(srt), " cells")
    } else {
      cm <- srt; message("[WARN] 无CM关键词，使用全部")
    }
  } else {
    cm <- srt; message("[WARN] 无细胞类型列，使用全部")
  }
}

# ========================================
# T26: NDUFAF3-NDUFB7 共表达
# ========================================
message("")
message("========================================")
message("T26: NDUFAF3-NDUFB7共表达验证")
message("========================================")

genes_t26 <- c("NDUFAF3", "NDUFB7")
avail_t26 <- intersect(genes_t26, rownames(cm))

if (length(avail_t26) < 2) {
  message("[FAIL] 缺失基因。可用: ", paste(avail_t26, collapse = ", "))
  all_g <- rownames(cm)
  for (g in genes_t26) {
    m <- grep(g, all_g, ignore.case = TRUE, value = TRUE)
    message("  匹配 '", g, "': ", paste(head(m, 5), collapse = ", "))
  }
  T26_PASS <- FALSE
} else {
  expr26 <- FetchData(cm, vars = avail_t26)
  colnames(expr26) <- avail_t26
  
  cp <- cor.test(expr26[[1]], expr26[[2]], method = "pearson")
  cs <- cor.test(expr26[[1]], expr26[[2]], method = "spearman")
  
  message("Pearson r  = ", round(cp$estimate, 3), "  p = ", format(cp$p.value, digits = 2, scientific = TRUE))
  message("Spearman ρ = ", round(cs$estimate, 3),     "  p = ", format(cs$p.value, digits = 2, scientific = TRUE))
  message("N = ", nrow(expr26))
  
  p26 <- ggplot(expr26, aes(x = .data[[avail_t26[1]]], y = .data[[avail_t26[2]]])) +
    geom_point(alpha = 0.3, size = 0.5, color = "#31688E") +
    geom_smooth(method = "lm", color = "#FDE725", se = TRUE, linewidth = 1) +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("r = ", round(cp$estimate, 3), "\np = ", format(cp$p.value, digits = 1, scientific = TRUE)),
             hjust = 1.1, vjust = -0.5, size = 3) +
    labs(title = "NDUFAF3-NDUFB7 Co-expression", subtitle = LABEL,
         x = avail_t26[1], y = avail_t26[2]) +
    theme_minimal(base_size = 10)
  
  outdir <- file.path(PROJECT_DIR, "03_results/T26_NDUFAF3_NDUFB7")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(outdir, "V92_T26_scatter.png"), p26, width = 6, height = 5, dpi = 300)
  
  stats26 <- data.frame(Data_Source = LABEL, N = nrow(expr26),
                        Pearson_r = cp$estimate, Pearson_p = cp$p.value,
                        Spearman_rho = cs$estimate, Spearman_p = cs$p.value)
  write.csv(stats26, file.path(outdir, "V92_T26_stats.csv"), row.names = FALSE)
  
  if (cp$estimate > 0.6 && cp$p.value < 0.05) {
    message("[PASS] 共表达显著且强相关！支持'组装缺陷轴'假说")
    T26_PASS <- TRUE
  } else if (cp$p.value < 0.05) {
    message("[PARTIAL] 共表达存在但强度中等")
    T26_PASS <- TRUE
  } else {
    message("[WARN] 共表达不显著")
    T26_PASS <- FALSE
  }
  message("[DONE] T26 -> ", outdir)
}

# ========================================
# T27: OXPHOS多复合体崩溃
# ========================================
message("")
message("========================================")
message("T27: OXPHOS多复合体协同崩溃验证")
message("========================================")

complex_genes <- list(
  Complex_I   = c("NDUFB7", "NDUFB8", "NDUFB10", "NDUFA9", "NDUFS1"),
  Complex_II  = c("SDHA", "SDHB", "SDHC", "SDHD"),
  Complex_III = c("UQCRC1", "UQCRC2", "CYTB", "CYC1"),
  Complex_IV  = c("COX4I1", "COX5A", "COX5B", "COX6C"),
  Complex_V   = c("ATP5F1A", "ATP5F1B", "ATP5MC1", "ATP5ME")
)

all_genes <- unlist(complex_genes)
avail_t27 <- intersect(all_genes, rownames(cm))
message("可用基因: ", length(avail_t27), "/", length(all_genes))
if (length(setdiff(all_genes, avail_t27)) > 0) {
  message("缺失: ", paste(setdiff(all_genes, avail_t27), collapse = ", "))
}

if (length(avail_t27) < 5) {
  message("[FAIL] 可用基因过少，T27阻塞")
  T27_PASS <- FALSE
} else {
  expr27 <- FetchData(cm, vars = avail_t27)
  
  grp_col <- intersect(c("condition","Condition","disease.state.ch1","group","cell_type"), colnames(cm@meta.data))[1]
  if (!is.na(grp_col)) {
    expr27$condition <- cm@meta.data[[grp_col]]
    message("[INFO] 分组 '", grp_col, "': ", paste(unique(expr27$condition), collapse = ", "))
  } else {
    expr27$condition <- "All"
    message("[WARN] 无分组列")
  }
  
  results <- data.frame()
  for (cname in names(complex_genes)) {
    g <- intersect(complex_genes[[cname]], colnames(expr27))
    if (length(g) == 0) next
    expr27[[cname]] <- rowMeans(expr27[, g, drop = FALSE])
    
    if (length(unique(expr27$condition)) >= 2) {
      grp <- split(expr27[[cname]], expr27$condition)
      if (length(grp) >= 2) {
        g1 <- grp[[1]]; g2 <- grp[[2]]
        tt <- tryCatch(t.test(g1, g2), error = function(e) NULL)
        if (!is.null(tt)) {
          results <- rbind(results, data.frame(
            Complex = cname, N_Genes = length(g),
            Group1 = names(grp)[1], Group1_Mean = mean(g1, na.rm = TRUE),
            Group2 = names(grp)[2], Group2_Mean = mean(g2, na.rm = TRUE),
            Log2FC = log2(mean(g2, na.rm = TRUE) / max(mean(g1, na.rm = TRUE), 1e-6)),
            P_Value = tt$p.value,
            Direction = ifelse(mean(g2, na.rm = TRUE) < mean(g1, na.rm = TRUE), "DOWN", "UP")
          ))
        }
      }
    }
  }
  
  outdir27 <- file.path(PROJECT_DIR, "03_results/T27_OXPHOS_Collapse")
  dir.create(outdir27, showWarnings = FALSE, recursive = TRUE)
  
  if (nrow(results) > 0) {
    write.csv(results, file.path(outdir27, "V92_T27_stats.csv"), row.names = FALSE)
    message(""); print(results[, c("Complex","Log2FC","P_Value","Direction")])
    
    down_n <- sum(results$Direction == "DOWN" & results$P_Value < 0.05, na.rm = TRUE)
    message(""); message("显著下调: ", down_n, "/", nrow(results))
    
    if (down_n >= 3) {
      message("[PASS] ≥3个复合体下调，支持'多复合体协同崩溃'")
      T27_PASS <- TRUE
    } else if (down_n >= 2) {
      message("[PARTIAL] 2个复合体下调")
      T27_PASS <- TRUE
    } else {
      message("[WARN] 下调不足，降级为'Complex I特异性'")
      T27_PASS <- FALSE
    }
  } else {
    message("[WARN] 无组间统计")
    T27_PASS <- FALSE
  }
  message("[DONE] T27 -> ", outdir27)
}

# ========================================
# 状态板
# ========================================
message("")
message("========================================")
message("V92 最终状态")
message("========================================")
message("数据源: ", LABEL)
t26_status <- ifelse(exists("T26_PASS") && T26_PASS, "✅ PASS", ifelse(exists("T26_PASS"), "❌ FAIL", "⏭️ SKIP"))
t27_status <- ifelse(exists("T27_PASS") && T27_PASS, "✅ PASS", ifelse(exists("T27_PASS"), "❌ FAIL", "⏭️ SKIP"))
message("T26 (NDUFAF3-NDUFB7): ", t26_status)
message("T27 (OXPHOS崩溃):     ", t27_status)
message("========================================")
