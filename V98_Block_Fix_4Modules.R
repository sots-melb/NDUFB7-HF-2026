#!/usr/bin/env Rscript
# V98: 四模块独立修复（R4→F4→R3→F3）
# 任一模块失败不阻塞其他，全部固定路径

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(reshape2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT_DIR)

message("========================================")
message("V98: BLOCKED项四模块修复")
message("========================================")

REPORTS <- list()
safe_read <- function(path, type = "csv") {
  fp <- file.path(PROJECT_DIR, path)
  if (!file.exists(fp)) return(NULL)
  if (type == "csv") return(read.csv(fp, stringsAsFactors = FALSE))
  if (type == "rds") return(readRDS(fp))
  NULL
}

# ========================================
# MODULE 1: R4 铁死亡独立性重建（最可行，49M数据）
# ========================================
message("")
message(">>> [MODULE 1] R4: 铁死亡/凋亡/焦亡独立性重建")

tryCatch({
  # 从V93已用过的CM文件重建（49M，加载快）
  CM_FILE <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"
  
  if (!file.exists(CM_FILE)) stop("CM文件不存在")
  
  cm <- readRDS(CM_FILE)
  message("  CM数据: ", ncol(cm), " cells")
  
  # 三套死亡基因集
  death_sets <- list(
    Ferroptosis = c("GPX4", "SLC7A11", "ACSL4", "LPCAT3", "FTH1", "FTL", "NFE2L2", "HMOX1"),
    Apoptosis   = c("BAX", "BAK1", "CASP3", "CASP8", "CASP9", "BCL2", "BCL2L1", "PARP1"),
    Pyroptosis  = c("NLRP3", "GSDMD", "GSDME", "IL1B", "IL18", "CASP1", "CASP4", "CASP5")
  )
  
  score_cols <- character(0)
  for (name in names(death_sets)) {
    genes <- intersect(death_sets[[name]], rownames(cm))
    if (length(genes) >= 3) {
      cm <- AddModuleScore(cm, features = list(genes), name = name)
      col_name <- paste0(name, "1")
      score_cols <- c(score_cols, col_name)
      message("  [SCORE] ", name, ": ", length(genes), " genes")
    } else {
      message("  [SKIP] ", name, ": 仅", length(genes), "个基因")
    }
  }
  
  if (length(score_cols) >= 2) {
    score_df <- cm@meta.data[, score_cols, drop = FALSE]
    # 重命名
    new_names <- gsub("1$", "", score_cols)
    colnames(score_df) <- new_names
    
    # 计算相关性
    cor_mat <- cor(score_df, use = "pairwise.complete.obs", method = "spearman")
    message("\n  === 死亡方式相关性矩阵 ===")
    print(round(cor_mat, 3))
    
    # 判断独立性
    ferro_apop <- ifelse("Ferroptosis" %in% rownames(cor_mat) && "Apoptosis" %in% colnames(cor_mat),
                         cor_mat["Ferroptosis", "Apoptosis"], NA)
    ferro_pyro <- ifelse("Ferroptosis" %in% rownames(cor_mat) && "Pyroptosis" %in% colnames(cor_mat),
                         cor_mat["Ferroptosis", "Pyroptosis"], NA)
    
    message("\n  Ferroptosis-Apoptosis ρ = ", round(ferro_apop, 3))
    message("  Ferroptosis-Pyroptosis ρ = ", round(ferro_pyro, 3))
    
    # 判决
    if (!is.na(ferro_apop) && !is.na(ferro_pyro)) {
      verdict <- ifelse(abs(ferro_apop) < 0.5 && abs(ferro_pyro) < 0.5, "INDEPENDENT",
                 ifelse(abs(ferro_apop) < 0.7, "PARTIALLY_INDEPENDENT", "OVERLAPPING"))
    } else if (!is.na(ferro_apop)) {
      verdict <- ifelse(abs(ferro_apop) < 0.5, "INDEPENDENT (仅Fer-Apop)", "OVERLAPPING")
    } else {
      verdict <- "INSUFFICIENT_DATA"
    }
    
    message("  [", verdict, "]")
    
    # 保存
    outdir <- file.path(PROJECT_DIR, "03_results/V98_R4_Death_Independence")
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    write.csv(score_df, file.path(outdir, "V98_death_scores.csv"), row.names = FALSE)
    write.csv(as.data.frame(cor_mat), file.path(outdir, "V98_death_correlation.csv"))
    
    # 可视化
    if (!is.na(ferro_apop)) {
      p <- ggplot(score_df, aes(x = Ferroptosis, y = Apoptosis)) +
        geom_point(alpha = 0.3, color = "#31688E", size = 0.5) +
        geom_smooth(method = "lm", color = "#FDE725", se = TRUE) +
        annotate("text", x = Inf, y = -Inf,
                 label = paste0("Spearman ρ = ", round(ferro_apop, 3)),
                 hjust = 1.1, vjust = -0.5, size = 3) +
        labs(title = "Ferroptosis vs Apoptosis Independence in CM",
             x = "Ferroptosis Module Score", y = "Apoptosis Module Score") +
        theme_minimal(base_size = 10)
      ggsave(file.path(outdir, "V98_ferro_vs_apop.png"), p, width = 6, height = 5, dpi = 300)
    }
    
    REPORTS[[length(REPORTS)+1]] <- list(
      Module = "R4", Verdict = verdict,
      Evidence = paste0("Fer-Apop=", round(ferro_apop, 2), " Fer-Pyro=", round(ferro_pyro, 2)),
      N = nrow(score_df)
    )
    message("  [DONE] R4修复完成 -> ", outdir)
  } else {
    message("  [FAIL] 可用评分不足")
    REPORTS[[length(REPORTS)+1]] <- list(Module = "R4", Verdict = "BLOCKED", Evidence = "基因不足", N = 0)
  }
  
}, error = function(e) {
  message("  [ERROR] R4: ", conditionMessage(e))
  REPORTS[[length(REPORTS)+1]] <- list(Module = "R4", Verdict = "ERROR", Evidence = conditionMessage(e), N = 0)
})

# ========================================
# MODULE 2: F4 NYHA剂量反应（检查V78/V71文件）
# ========================================
message("")
message(">>> [MODULE 2] F4: NYHA剂量反应修复")

tryCatch({
  # 尝试多个候选文件（固定路径，不搜索）
  candidates <- c(
    "V78_GSE135055_NDUFB7_clinical.csv",
    "V71_GSE135055_full_phenotype_SOFT.csv",
    "V71_GSE135055_RNAseq_phenotype.csv",
    "V70_GSE135055_phenotype.csv"
  )
  
  best <- NULL
  for (c in candidates) {
    if (file.exists(c)) {
      tmp <- read.csv(c, stringsAsFactors = FALSE)
      message("  检查: ", basename(c), " (", ncol(tmp), " cols, ", nrow(tmp), " rows)")
      
      # 检查是否有NYHA和NDUFB7
      has_nyha <- any(grepl("NYHA|nyha|class", colnames(tmp), ignore.case = TRUE))
      has_ndufb7 <- any(grepl("NDUFB7|ndufb7", colnames(tmp), ignore.case = TRUE))
      
      if (has_nyha && has_ndufb7) {
        best <- tmp
        message("  [FOUND] ", basename(c), " 含NYHA+NDUFB7")
        break
      } else if (has_nyha) {
        message("  [PARTIAL] ", basename(c), " 有NYHA但无NDUFB7")
      }
    }
  }
  
  if (is.null(best)) {
    # 尝试从表达矩阵+SOFT表型手动合并
    pheno <- safe_read("V71_GSE135055_full_phenotype_SOFT.csv")
    expr <- safe_read("GSE135055_FPKM_Expression_Matrix.txt.gz")
    
    if (!is.null(pheno) && !is.null(expr)) {
      message("  [INFO] 尝试从表达矩阵提取NDUFB7并关联SOFT表型...")
      # 简化：假设表达矩阵第一列是基因名
      gene_col <- colnames(expr)[1]
      ndufb7_row <- grep("NDUFB7", expr[[gene_col]], ignore.case = TRUE)
      
      if (length(ndufb7_row) > 0) {
        ndufb7_vals <- as.numeric(expr[ndufb7_row[1], -1])
        sample_ids <- colnames(expr)[-1]
        
        # 找pheno中的样本ID列和NYHA列
        id_col <- grep("sample|gsm|title", colnames(pheno), ignore.case = TRUE)[1]
        nyha_col <- grep("NYHA|nyha|class|stage", colnames(pheno), ignore.case = TRUE)[1]
        
        if (!is.na(id_col) && !is.na(nyha_col)) {
          # 匹配（可能需要模糊匹配）
          pheno_ids <- as.character(pheno[[id_col]])
          matched <- sapply(sample_ids, function(s) {
            idx <- grep(s, pheno_ids, ignore.case = TRUE)[1]
            ifelse(is.na(idx), NA, idx)
          })
          
          valid <- !is.na(matched)
          if (sum(valid) >= 10) {
            nyha_vals <- pheno[[nyha_col]][matched[valid]]
            ndufb7_matched <- ndufb7_vals[valid]
            
            # 清理NYHA
            nyha_clean <- gsub("NYHA class |NYHA |Class ", "", nyha_vals, ignore.case = TRUE)
            nyha_num <- suppressWarnings(as.numeric(nyha_clean))
            
            valid2 <- !is.na(nyha_num) & nyha_num >= 1 & nyha_num <= 4 & !is.na(ndufb7_matched)
            if (sum(valid2) >= 10) {
              nyha_final <- nyha_num[valid2]
              ndufb7_final <- ndufb7_matched[valid2]
              
              sp <- cor.test(ndufb7_final, nyha_final, method = "spearman")
              kt <- kruskal.test(ndufb7_final ~ factor(nyha_final))
              
              verdict <- ifelse(sp$estimate < -0.2 && sp$p.value < 0.05, "STRONG_SUPPORT",
                         ifelse(sp$p.value < 0.05, "MODERATE_SUPPORT", "NO_SUPPORT"))
              
              message("\n  [", verdict, "] NYHA vs NDUFB7:")
              message("    Spearman ρ = ", round(sp$estimate, 3), " p = ", format(sp$p.value, digits = 2))
              message("    Kruskal-Wallis p = ", format(kt$p.value, digits = 2))
              message("    N = ", sum(valid2))
              
              outdir <- file.path(PROJECT_DIR, "03_results/V98_F4_NYHA")
              dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
              
              res_df <- data.frame(
                NYHA = nyha_final, NDUFB7 = ndufb7_final,
                Sample = sample_ids[valid][valid2]
              )
              write.csv(res_df, file.path(outdir, "V98_NYHA_NDUFB7.csv"), row.names = FALSE)
              
              # 箱线图
              p <- ggplot(res_df, aes(x = factor(NYHA), y = NDUFB7, fill = factor(NYHA))) +
                geom_boxplot(alpha = 0.7) +
                geom_jitter(alpha = 0.3, width = 0.2, size = 0.5) +
                labs(title = "NDUFB7 Expression by NYHA Class",
                     subtitle = paste0("Spearman ρ=", round(sp$estimate, 3), 
                                      ", p=", format(sp$p.value, digits = 2)),
                     x = "NYHA Class", y = "NDUFB7 Expression") +
                theme_minimal(base_size = 10) +
                theme(legend.position = "none")
              ggsave(file.path(outdir, "V98_NYHA_boxplot.png"), p, width = 6, height = 4, dpi = 300)
              
              REPORTS[[length(REPORTS)+1]] <- list(
                Module = "F4", Verdict = verdict,
                Evidence = paste0("ρ=", round(sp$estimate, 2), " p=", format(sp$p.value, digits = 2)),
                N = sum(valid2)
              )
              message("  [DONE] F4修复完成 -> ", outdir)
            } else {
              message("  [FAIL] 有效样本不足或NYHA解析失败")
              REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "BLOCKED", Evidence = "NYHA解析失败", N = 0)
            }
          } else {
            message("  [FAIL] 样本匹配不足")
            REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "BLOCKED", Evidence = "样本匹配<10", N = 0)
          }
        } else {
          message("  [FAIL] 找不到样本ID或NYHA列")
          REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "BLOCKED", Evidence = "列名不匹配", N = 0)
        }
      } else {
        message("  [FAIL] NDUFB7不在表达矩阵中")
        REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "BLOCKED", Evidence = "NDUFB7缺失", N = 0)
      }
    } else {
      message("  [FAIL] 表型或表达矩阵缺失")
      REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "BLOCKED", Evidence = "文件缺失", N = 0)
    }
  }
  
}, error = function(e) {
  message("  [ERROR] F4: ", conditionMessage(e))
  REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "ERROR", Evidence = conditionMessage(e), N = 0)
})

# ========================================
# MODULE 3: R3 LVAD恢复（尝试6.4G文件）
# ========================================
message("")
message(">>> [MODULE 3] R3: LVAD治疗后恢复验证")

tryCatch({
  GLOBAL_FILE <- "01_data/01_raw_geo/GSE226314_global.rds.gz"
  
  if (!file.exists(GLOBAL_FILE)) {
    # 检查Downloads
    GLOBAL_FILE <- "Downloads/GSE226314_global.rds.gz"
  }
  
  if (!file.exists(GLOBAL_FILE)) {
    message("  [MISS] GSE226314_global.rds.gz不存在")
    REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "文件缺失", N = 0)
  } else {
    sz <- file.info(GLOBAL_FILE)$size
    message("  [FOUND] ", basename(GLOBAL_FILE), " (", round(sz/1e9, 1), " GB)")
    message("  [LOAD] 尝试读取（可能需要2-5分钟，内存>8G）...")
    
    # 尝试读取
    con <- gzfile(GLOBAL_FILE, "rb")
    obj <- readRDS(con)
    close(con)
    
    message("  [PASS] 加载成功，对象类型: ", paste(class(obj), collapse = ", "))
    
    # 提取NDUFB7和pre/post分组
    if ("Seurat" %in% class(obj)) {
      # 找分组列
      meta <- obj@meta.data
      grp_col <- intersect(c("condition", "Condition", "time_point", "treatment", "LVAD", "group"), colnames(meta))[1]
      
      if (!is.na(grp_col)) {
        grps <- unique(meta[[grp_col]])
        message("  分组列 '", grp_col, "': ", paste(grps, collapse = ", "))
        
        # 找pre/post关键词
        pre_mask <- grepl("pre|baseline|before|0h", meta[[grp_col]], ignore.case = TRUE)
        post_mask <- grepl("post|after|LVAD|1w|4w", meta[[grp_col]], ignore.case = TRUE)
        
        if (sum(pre_mask) > 0 && sum(post_mask) > 0) {
          # 提取NDUFB7
          if ("NDUFB7" %in% rownames(obj)) {
            expr <- FetchData(obj, vars = "NDUFB7")
            pre_vals <- expr$NDUFB7[pre_mask]
            post_vals <- expr$NDUFB7[post_mask]
            
            tt <- t.test(pre_vals, post_vals)
            delta <- mean(post_vals) - mean(pre_vals)
            pct <- delta / mean(pre_vals) * 100
            
            verdict <- ifelse(delta > 0 && tt$p.value < 0.05, "REVERSED",
                       ifelse(delta > 0, "PARTIALLY_REVERSED", "NOT_REVERSED"))
            
            message("\n  [", verdict, "]")
            message("    Pre-LVAD  mean: ", round(mean(pre_vals), 4))
            message("    Post-LVAD mean: ", round(mean(post_vals), 4))
            message("    Δ = ", round(delta, 4), " (", round(pct, 1), "%)")
            message("    t-test p = ", format(tt$p.value, digits = 2))
            message("    N_pre = ", length(pre_vals), ", N_post = ", length(post_vals))
            
            outdir <- file.path(PROJECT_DIR, "03_results/V98_R3_LVAD")
            dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
            
            res_df <- data.frame(
              Group = c(rep("Pre-LVAD", length(pre_vals)), rep("Post-LVAD", length(post_vals))),
              NDUFB7 = c(pre_vals, post_vals)
            )
            write.csv(res_df, file.path(outdir, "V98_LVAD_NDUFB7.csv"), row.names = FALSE)
            
            p <- ggplot(res_df, aes(x = Group, y = NDUFB7, fill = Group)) +
              geom_boxplot(alpha = 0.7) +
              geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
              annotate("text", x = 1.5, y = Inf, 
                       label = paste0("Δ=", round(pct, 1), "%, p=", format(tt$p.value, digits = 1)),
                       vjust = 1.5, size = 3) +
              labs(title = "NDUFB7 Recovery After LVAD",
                   subtitle = paste0("N_pre=", length(pre_vals), ", N_post=", length(post_vals))) +
              theme_minimal(base_size = 10) +
              theme(legend.position = "none")
            ggsave(file.path(outdir, "V98_LVAD_boxplot.png"), p, width = 5, height = 4, dpi = 300)
            
            REPORTS[[length(REPORTS)+1]] <- list(
              Module = "R3", Verdict = verdict,
              Evidence = paste0("Δ=", round(pct, 1), "% p=", format(tt$p.value, digits = 2)),
              N = paste0(length(pre_vals), "/", length(post_vals))
            )
            message("  [DONE] R3修复完成 -> ", outdir)
            
          } else {
            message("  [FAIL] NDUFB7不在基因列表中")
            REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "NDUFB7缺失", N = 0)
          }
        } else {
          message("  [FAIL] 未识别pre/post分组")
          message("  可用分组值: ", paste(grps, collapse = ", "))
          REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "分组不匹配", N = 0)
        }
      } else {
        message("  [FAIL] 无分组列")
        message("  可用列: ", paste(head(colnames(meta), 10), collapse = ", "))
        REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "无分组列", N = 0)
      }
    } else {
      message("  [FAIL] 非Seurat对象")
      REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "类型错误", N = 0)
    }
  }
  
}, error = function(e) {
  msg <- conditionMessage(e)
  message("  [ERROR] R3: ", msg)
  if (grepl("cannot allocate|memory|vector", msg, ignore.case = TRUE)) {
    message("  [DIAGNOSIS] 内存不足，无法加载6.4G文件")
    message("  [WORKAROUND] 建议:")
    message("    1. 关闭其他R进程释放内存")
    message("    2. 在服务器上执行（非本地PC）")
    message("    3. 或在Discussion中写'Future work: validate in independent LVAD cohort'")
  }
  REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "ERROR", Evidence = msg, N = 0)
})

# ========================================
# MODULE 4: F3 GSE243655 Fer-1解压+解析
# ========================================
message("")
message(">>> [MODULE 4] F3: GSE243655 Fer-1数据准备")

tryCatch({
  RAW_TAR <- "01_data/01_raw_geo/GSE243655/GSE243655_RAW(1).tar"
  if (!file.exists(RAW_TAR)) {
    RAW_TAR <- "Downloads/GSE243655_RAW(1).tar"
  }
  
  if (!file.exists(RAW_TAR)) {
    message("  [MISS] GSE243655_RAW.tar不存在")
    REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "BLOCKED", Evidence = "文件缺失", N = 0)
  } else {
    message("  [FOUND] ", basename(RAW_TAR))
    
    # 解压目录
    EXTRACT_DIR <- file.path(PROJECT_DIR, "01_data/01_raw_geo/GSE243655")
    dir.create(EXTRACT_DIR, showWarnings = FALSE, recursive = TRUE)
    
    # 检查是否已解压
    existing <- list.files(EXTRACT_DIR, pattern = "\\.txt$|\\.csv$|\\.cel$", ignore.case = TRUE)
    
    if (length(existing) > 0) {
      message("  [INFO] 已解压文件: ", paste(existing, collapse = ", "))
    } else {
      message("  [EXTRACT] 解压中...")
      system(paste("tar -xf", shQuote(RAW_TAR), "-C", shQuote(EXTRACT_DIR)), intern = TRUE)
      existing <- list.files(EXTRACT_DIR, pattern = "\\.txt$|\\.csv$|\\.cel$", ignore.case = TRUE)
      message("  [DONE] 解压完成，文件数: ", length(existing))
    }
    
    if (length(existing) > 0) {
      message("  [INFO] 文件列表:")
      for (f in head(existing, 5)) {
        message("    - ", f)
      }
      
      # 尝试读取series_matrix获取分组
      sm <- safe_read("01_data/01_raw_geo/GSE243655/GSE243655_series_matrix.txt.gz")
      if (is.null(sm)) sm <- safe_read("Downloads/GSE243655_series_matrix.txt.gz")
      
      if (!is.null(sm)) {
        message("  [INFO] Series matrix已加载，需手动检查Fer-1 vs Vehicle分组")
        message("  [PENDING] F3数据已就绪，需手动:")
        message("    1. 确认GSE243655分组（Fer-1 vs Vehicle）")
        message("    2. 提取NDUFB7表达值")
        message("    3. 做t-test比较")
      } else {
        message("  [PENDING] 数据已解压，但series matrix未找到")
      }
      
      REPORTS[[length(REPORTS)+1]] <- list(
        Module = "F3", Verdict = "PENDING_MANUAL",
        Evidence = paste0("解压完成，", length(existing), "个文件待解析"),
        N = length(existing)
      )
      message("  [DONE] F3数据准备完成，等待手动解析分组")
    } else {
      message("  [FAIL] 解压后无识别文件")
      REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "BLOCKED", Evidence = "解压失败", N = 0)
    }
  }
  
}, error = function(e) {
  message("  [ERROR] F3: ", conditionMessage(e))
  REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "ERROR", Evidence = conditionMessage(e), N = 0)
})

# ========================================
# 综合报告
# ========================================
message("")
message("========================================")
message("V98 修复结果汇总")
message("========================================")

report_df <- do.call(rbind, lapply(REPORTS, function(x) {
  data.frame(Module = x$Module, Verdict = x$Verdict, Evidence = x$Evidence, N = x$N, stringsAsFactors = FALSE)
}))
print(report_df, row.names = FALSE)

outdir <- file.path(PROJECT_DIR, "03_results/V98_Validation_Fix")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(report_df, file.path(outdir, "V98_fix_report.csv"), row.names = FALSE)

message("")
message("[DONE] V98修复完成。查看各模块结果:")
message("  R4: 03_results/V98_R4_Death_Independence/")
message("  F4: 03_results/V98_F4_NYHA/")
message("  R3: 03_results/V98_R3_LVAD/")
message("  F3: 01_data/01_raw_geo/GSE243655/ (待手动分组)")
