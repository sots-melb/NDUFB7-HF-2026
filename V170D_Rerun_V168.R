suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ppcor)
})
PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V170D: V168B Complex I + V168D Partial Cor 重跑")
message("========================================")

# 加载抢救的矩阵
rescue_file <- "03_results/V170_Rescue/GSE57338_rescued_matrix.rds"
if(!file.exists(rescue_file)) {
  message("[FAIL] Run V170A first")
  quit(save="no", status=1)
}

expr_mat <- readRDS(rescue_file)
message("[1/4] Loaded: ", nrow(expr_mat), " x ", ncol(expr_mat))

# 确保行名可用
if(is.null(rownames(expr_mat))) {
  message("[FAIL] Matrix lacks rownames")
  quit(save="no", status=1)
}

outdir <- "03_results/V170_Rescue/V168_Rerun"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ========================================
# V168B: Complex I 特异性
# ========================================
message("\n[2/4] V168B: Complex I vs Pan-Death correlation...")

c1_genes <- c("NDUFB7", "NDUFA13", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB8", 
              "NDUFB9", "NDUFB10", "NDUFB11", "NDUFA1", "NDUFA2", "NDUFA3",
              "NDUFS1", "NDUFS2", "NDUFS3", "NDUFV1", "NDUFV2")

found_c1 <- c1_genes[tolower(c1_genes) %in% tolower(rownames(expr_mat))]
message("  Complex I found: ", length(found_c1), "/", length(c1_genes))
if(length(found_c1) < 3) {
  message("[SKIP] Too few Complex I genes")
} else {
  # 死亡通路评分
  pw_defs <- list(
    Apoptosis = c("BCL2", "BAX", "CASP3", "CASP7", "CASP9", "BAD", "BAK1"),
    Necroptosis = c("RIPK1", "RIPK3", "MLKL", "CYLD", "TNF"),
    Ferroptosis_Defense = c("GPX4", "SLC7A11", "FTH1", "FTL", "NFE2L2"),
    Ferroptosis_Execution = c("ACSL4", "LPCAT3", "LOX", "PTGS2"),
    Autophagy = c("ATG5", "ATG7", "BECN1", "MAP1LC3B", "SQSTM1"),
    Pyroptosis = c("NLRP3", "GSDMD", "IL1B", "IL18", "CASP1")
  )
  
  pw_scores <- data.frame(Sample = colnames(expr_mat), stringsAsFactors = FALSE)
  for(pw in names(pw_defs)) {
    g <- pw_defs[[pw]][tolower(pw_defs[[pw]]) %in% tolower(rownames(expr_mat))]
    if(length(g) > 0) {
      vals <- expr_mat[tolower(g), , drop = FALSE]
      if(is.matrix(vals) && nrow(vals) > 0) {
        pw_scores[[pw]] <- colMeans(vals, na.rm = TRUE)
        message("  ", pw, ": ", length(g), " genes")
      }
    }
  }
  
  # 计算Complex I与通路的相关
  c1_mat <- expr_mat[tolower(found_c1), , drop = FALSE]
  results <- data.frame()
  for(g in rownames(c1_mat)) {
    for(j in 2:ncol(pw_scores)) {
      pw <- names(pw_scores)[j]
      r <- cor(c1_mat[g, ], pw_scores[[j]], method = "spearman", use = "pairwise.complete.obs")
      if(is.finite(r)) {
        results <- rbind(results, data.frame(Gene = g, Pathway = pw, Rho = round(r, 3), stringsAsFactors = FALSE))
      }
    }
  }
  
  # NDUFB7特异性判定
  nduf_res <- results[results$Gene == "NDUFB7", ]
  others_res <- results[results$Gene != "NDUFB7", ]
  
  summary_b <- data.frame()
  for(pw in unique(results$Pathway)) {
    nduf_rho <- nduf_res$Rho[nduf_res$Pathway == pw]
    other_rho <- others_res$Rho[others_res$Pathway == pw]
    if(length(nduf_rho) > 0 && length(other_rho) > 0) {
      summary_b <- rbind(summary_b, data.frame(
        Pathway = pw,
        NDUFB7_Rho = round(nduf_rho, 3),
        Max_Other_Rho = round(max(abs(other_rho), na.rm = TRUE), 3),
        Median_Other_Rho = round(median(abs(other_rho), na.rm = TRUE), 3),
        NDUFB7_Is_Max = abs(nduf_rho) >= max(abs(other_rho), na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  fwrite(summary_b, file.path(outdir, "V168B_Complex1_Specificity.csv"))
  message("\n=== Complex I Specificity ===")
  print(summary_b)
  
  n_max <- sum(summary_b$NDUFB7_Is_Max, na.rm = TRUE)
  message("\n[RESULT] NDUFB7 is STRONGEST in ", n_max, "/", nrow(summary_b), " pathways")
  if(n_max >= 4) {
    message("[PASS] SPECIFIC pan-death signature confirmed")
  } else if(n_max >= 2) {
    message("[PARTIAL] Among top correlates")
  } else {
    message("[WARN] Not specific — other Complex I genes stronger")
  }
}

# ========================================
# V168D: 偏相关（控制ROS）
# ========================================
message("\n[3/4] V168D: Partial correlation controlling ROS...")

ros_genes <- c("NOX4", "SOD2", "CAT", "GPX1")
death_defs <- list(
  Apoptosis = c("BCL2", "BAX", "CASP3", "CASP7", "CASP9"),
  Necroptosis = c("RIPK1", "RIPK3", "MLKL"),
  Ferroptosis_Defense = c("GPX4", "SLC7A11", "FTH1", "FTL"),
  Ferroptosis_Execution = c("ACSL4", "LPCAT3", "LOX")
)

all_needed <- unique(c("NDUFB7", ros_genes, unlist(death_defs)))
found_all <- all_needed[tolower(all_needed) %in% tolower(rownames(expr_mat))]
message("  Target genes found: ", length(found_all), "/", length(all_needed))

ros_found <- ros_genes[ros_genes %in% found_all]
if(length(ros_found) < 2) {
  message("[SKIP] Too few ROS genes (", paste(ros_found, collapse = ", "), ")")
} else {
  sub_mat <- expr_mat[tolower(found_all), , drop = FALSE]
  
  # 计算样本级通路评分
  pw_scores2 <- data.frame(Sample = colnames(sub_mat), stringsAsFactors = FALSE)
  for(pw in names(death_defs)) {
    g <- death_defs[[pw]][tolower(death_defs[[pw]]) %in% tolower(rownames(sub_mat))]
    if(length(g) > 0) {
      vals <- sub_mat[tolower(g), , drop = FALSE]
      pw_scores2[[pw]] <- colMeans(vals, na.rm = TRUE)
    }
  }
  
  # 构建数据框
  cor_df <- data.frame(NDUFB7 = as.numeric(sub_mat["NDUFB7", ]), stringsAsFactors = FALSE)
  for(pw in setdiff(names(pw_scores2), "Sample")) {
    cor_df[[pw]] <- as.numeric(pw_scores2[[pw]])
  }
  for(rg in ros_found) {
    cor_df[[rg]] <- as.numeric(sub_mat[rg, ])
  }
  
  cor_df <- cor_df[complete.cases(cor_df), ]
  message("  Complete cases: ", nrow(cor_df), " samples")
  
  if(nrow(cor_df) > 50) {
    results_pc <- data.frame()
    for(pw in setdiff(names(pw_scores2), "Sample")) {
      if(pw %in% names(cor_df)) {
        tryCatch({
          pc <- pcor.test(cor_df$NDUFB7, cor_df[[pw]], cor_df[, ros_found, drop = FALSE])
          results_pc <- rbind(results_pc, data.frame(
            Pathway = pw,
            ROS_Control = paste(ros_found, collapse = "+"),
            Partial_R = round(pc$estimate, 3),
            Partial_P = format(pc$p.value, digits = 2, scientific = TRUE),
            N = nrow(cor_df),
            stringsAsFactors = FALSE
          ))
        }, error = function(e) message("  ", pw, " failed: ", conditionMessage(e)))
      }
    }
    
    if(nrow(results_pc) > 0) {
      # 与简单相关比较
      simple <- data.frame()
      for(pw in setdiff(names(pw_scores2), "Sample")) {
        if(pw %in% names(cor_df)) {
          r <- cor(cor_df$NDUFB7, cor_df[[pw]], method = "spearman")
          simple <- rbind(simple, data.frame(Pathway = pw, Simple_Rho = round(r, 3), stringsAsFactors = FALSE))
        }
      }
      
      comp <- merge(results_pc[, c("Pathway", "Partial_R")], simple, by = "Pathway")
      comp$Attenuation_Pct <- round((1 - abs(comp$Partial_R) / abs(comp$Simple_Rho)) * 100, 1)
      
      fwrite(comp, file.path(outdir, "V168D_Partial_Correlation.csv"))
      message("\n=== Partial Correlation (ROS-controlled) ===")
      print(comp[order(comp$Attenuation_Pct, decreasing = TRUE), ])
      
      high_att <- comp$Attenuation_Pct[comp$Attenuation_Pct > 50 & !is.na(comp$Attenuation_Pct)]
      if(length(high_att) > 0) {
        message("\n[WARN] ", length(high_att), " pathways show >50% attenuation")
        message("  → NDUFB7 may mediate death primarily through ROS")
      } else {
        message("\n[PASS] NDUFB7-death associations largely ROS-independent")
      }
    }
  }
}

# ========================================
# V168C: GSE59867 ACSL4/GPX4（如果V170B成功）
# ========================================
message("\n[4/4] V168C: GSE59867 ACSL4/GPX4 ratio validation...")
mg_file <- "03_results/V170_Rescue/GSE59867_multigene_extracted.csv"
if(file.exists(mg_file)) {
  mg <- fread(mg_file)
  message("  Loaded multigene data: ", nrow(mg), " genes")
  
  # 判断格式
  gene_col <- names(mg)[1]
  if(gene_col %in% c("SYMBOL", "PROBEID", "V1")) gene_col <- names(mg)[2]
  
  # 找ACSL4和GPX4行
  gene_names <- as.character(mg[[names(mg)[1]]])
  acsl4_idx <- which(tolower(gene_names) == "acsl4")
  gpx4_idx <- which(tolower(gene_names) == "gpx4")
  
  if(length(acsl4_idx) > 0 && length(gpx4_idx) > 0) {
    sample_cols <- setdiff(names(mg), c("PROBEID", "SYMBOL", "GENENAME", "ID_REF", "V1"))
    acsl4_vals <- as.numeric(mg[acsl4_idx, ..sample_cols])
    gpx4_vals <- as.numeric(mg[gpx4_idx, ..sample_cols])
    
    ratio <- log2((acsl4_vals + 1) / (gpx4_vals + 1))
    message("  ACSL4/GPX4 log2 ratio: ", round(mean(ratio, na.rm = TRUE), 3), 
            " ± ", round(sd(ratio, na.rm = TRUE), 3))
    
    # 保存
    ratio_df <- data.frame(Sample = sample_cols, ACSL4_GPX4_log2Ratio = ratio, stringsAsFactors = FALSE)
    fwrite(ratio_df, file.path(outdir, "V168C_ACSL4_GPX4_ratio.csv"))
    message("  [PASS] Saved ratio data")
  } else {
    message("  [SKIP] ACSL4 or GPX4 not found in multigene extract")
  }
} else {
  message("  [SKIP] V170B multigene file not available")
}

message("\n[DONE] V170D: All outputs in ", outdir)
