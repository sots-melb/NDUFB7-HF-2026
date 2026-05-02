#!/usr/bin/env Rscript
# V89: 加载正确的GSE183852文件 + T26共表达 + T27 OXPHOS崩溃

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(reshape2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

message("========================================")
message("V89: 正确文件加载 + T26/T27 连跑")
message("========================================")

# --- 策略：优先用 Pure_CM_Subsampled (915MB, 已subset) ---
# 备选：GSE183852_DCM_Nuclei.Robj (8.3GB, 全量)

CANDIDATES <- c(
  "Downloads/GSE183852_Pure_CM_Subsampled.RDS",
  "01_data/01_raw_geo/GSE183852_Pure_CM_Subsampled.RDS",
  "Downloads/GSE183852_DCM_Nuclei.Robj",
  "Downloads/GSE183852_DCM_Nuclei.Robj.gz",
  "01_data/01_raw_geo/GSE183852_DCM_Nuclei.Robj"
)

FILE <- NULL
for (c in CANDIDATES) {
  if (file.exists(c)) {
    FILE <- c
    message("[FOUND] 数据文件: ", basename(c), " (", round(file.info(c)$size/1024/1024, 1), " MB)")
    break
  }
}

if (is.null(FILE)) stop("[FAIL] 未找到任何GSE183852数据文件。请确认下载完成。")

# --- 加载 ---
message("[LOAD] 读取中，请等待...")
obj <- readRDS(FILE)
cls <- class(obj)

if ("Seurat" %in% cls) {
  srt <- obj
  message("[PASS] 标准Seurat对象，", ncol(srt), " cells")
} else if (is.list(obj) && any(sapply(obj, function(x) "Seurat" %in% class(x)))) {
  idx <- which(sapply(obj, function(x) "Seurat" %in% class(x)))[1]
  srt <- obj[[idx]]
  message("[PASS] List中提取Seurat，", ncol(srt), " cells")
} else {
  stop("[FAIL] 无法识别Seurat结构: ", paste(cls, collapse = ", "))
}

# --- 确认CM子集 ---
if (grepl("Pure_CM|CM", FILE, ignore.case = TRUE)) {
  cm <- srt
  message("[INFO] 文件已是Pure_CM，直接使用")
} else {
  # 从全量中提取CM
  meta_cols <- colnames(cm@meta.data)
  cm_col <- intersect(c("cell_type", "predicted.id", "celltype"), meta_cols)[1]
  if (!is.na(cm_col)) {
    cm_mask <- grepl("Cardio|CM|myocyte", srt@meta.data[[cm_col]], ignore.case = TRUE)
    if (sum(cm_mask) > 0) {
      cm <- subset(srt, cells = colnames(srt)[cm_mask])
      message("[SUBSET] CM子集: ", ncol(cm), " / ", ncol(srt), " cells")
    } else {
      cm <- srt
      message("[WARN] 未匹配CM，使用全部")
    }
  } else {
    cm <- srt
    message("[WARN] 无细胞类型列，使用全部")
  }
}

# ========================================
# T26: NDUFAF3-NDUFB7共表达验证
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
  
  # 可视化
  p26 <- ggplot(expr26, aes(x = .data[[avail_t26[1]]], y = .data[[avail_t26[2]]])) +
    geom_point(alpha = 0.3, size = 0.5, color = "#31688E") +
    geom_smooth(method = "lm", color = "#FDE725", se = TRUE, linewidth = 1) +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("r = ", round(cp$estimate, 3), "\np = ", format(cp$p.value, digits = 1, scientific = TRUE)),
             hjust = 1.1, vjust = -0.5, size = 3) +
    labs(title = "NDUFAF3-NDUFB7 Co-expression in CM",
         x = avail_t26[1], y = avail_t26[2]) +
    theme_minimal(base_size = 10)
  
  outdir <- file.path(PROJECT_DIR, "03_results/T26_NDUFAF3_NDUFB7")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(outdir, "V89_T26_scatter.png"), p26, width = 6, height = 5, dpi = 300)
  
  stats26 <- data.frame(
    N = nrow(expr26), Pearson_r = cp$estimate, Pearson_p = cp$p.value,
    Spearman_rho = cs$estimate, Spearman_p = cs$p.value
  )
  write.csv(stats26, file.path(outdir, "V89_T26_stats.csv"), row.names = FALSE)
  
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
  message("[DONE] T26结果: ", outdir)
}

# ========================================
# T27: OXPHOS多复合体崩溃验证
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
message("缺失: ", paste(setdiff(all_genes, avail_t27), collapse = ", "))

if (length(avail_t27) < 5) {
  message("[FAIL] 可用基因过少，T27阻塞")
  T27_PASS <- FALSE
} else {
  expr27 <- FetchData(cm, vars = avail_t27)
  
  # 获取分组（自适应）
  grp_col <- intersect(c("condition", "Condition", "disease.state.ch1", "disease_state", "group"), colnames(cm@meta.data))[1]
  
  if (!is.na(grp_col)) {
    expr27$condition <- cm@meta.data[[grp_col]]
    message("[INFO] 分组列: '", grp_col, "' | 组: ", paste(unique(expr27$condition), collapse = ", "))
  } else {
    expr27$condition <- "All"
    message("[WARN] 无分组列，仅输出描述统计")
  }
  
  # 逐复合体统计
  results <- data.frame()
  for (cname in names(complex_genes)) {
    g <- intersect(complex_genes[[cname]], colnames(expr27))
    if (length(g) == 0) next
    expr27[[cname]] <- rowMeans(expr27[, g, drop = FALSE])
    
    if (length(unique(expr27$condition)) >= 2) {
      grp <- split(expr27[[cname]], expr27$condition)
      if (length(grp) >= 2) {
        g1 <- grp[[1]]; g2 <- grp[[2]]
        tt <- t.test(g1, g2)
        results <- rbind(results, data.frame(
          Complex = cname, N_Genes = length(g),
          Group1 = names(grp)[1], Group1_Mean = mean(g1, na.rm = TRUE),
          Group2 = names(grp)[2], Group2_Mean = mean(g2, na.rm = TRUE),
          Log2FC = log2(mean(g2, na.rm = TRUE) / mean(g1, na.rm = TRUE)),
          P_Value = tt$p.value,
          Direction = ifelse(mean(g2, na.rm = TRUE) < mean(g1, na.rm = TRUE), "DOWN", "UP")
        ))
      }
    }
  }
  
  outdir27 <- file.path(PROJECT_DIR, "03_results/T27_OXPHOS_Collapse")
  dir.create(outdir27, showWarnings = FALSE, recursive = TRUE)
  
  if (nrow(results) > 0) {
    write.csv(results, file.path(outdir27, "V89_T27_stats.csv"), row.names = FALSE)
    message("")
    print(results[, c("Complex", "Log2FC", "P_Value", "Direction")])
    
    down_n <- sum(results$Direction == "DOWN" & results$P_Value < 0.05)
    message("")
    message("显著下调复合体: ", down_n, "/", nrow(results))
    
    if (down_n >= 3) {
      message("[PASS] ≥3个复合体显著下调，支持'多复合体协同崩溃'假说")
      T27_PASS <- TRUE
    } else if (down_n >= 2) {
      message("[PARTIAL] 2个复合体下调")
      T27_PASS <- TRUE
    } else {
      message("[WARN] 下调不足，叙事需降级为'Complex I特异性'")
      T27_PASS <- FALSE
    }
    
    # 热图
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      library(pheatmap)
      hm_mat <- as.matrix(results[, c("Group1_Mean", "Group2_Mean")])
      rownames(hm_mat) <- results$Complex
      colnames(hm_mat) <- c(results$Group1[1], results$Group2[1])
      
      png(file.path(outdir27, "V89_T27_heatmap.png"), width = 600, height = 400)
      pheatmap(hm_mat, scale = "row", main = "OXPHOS Complex Expression (z-score)",
               color = colorRampPalette(c("#440154", "#31688E", "#35B779", "#FDE725"))(100))
      dev.off()
      message("[DONE] 热图保存")
    }
    
  } else {
    message("[WARN] 无组间比较结果")
    T27_PASS <- FALSE
  }
  message("[DONE] T27结果: ", outdir27)
}

# ========================================
# 最终状态板
# ========================================
message("")
message("========================================")
message("V89 最终状态")
message("========================================")
message("T26 (NDUFAF3-NDUFB7): ", ifelse(T26_PASS, "✅ PASS", ifelse(exists("T26_PASS"), "❌ FAIL", "⏭️ SKIP")))
message("T27 (OXPHOS崩溃):     ", ifelse(T27_PASS, "✅ PASS", ifelse(exists("T27_PASS"), "❌ FAIL", "⏭️ SKIP")))
message("========================================")
