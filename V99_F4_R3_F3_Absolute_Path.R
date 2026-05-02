#!/usr/bin/env Rscript
# V99: F4/R3/F3 绝对路径修复
# 使用~/Downloads/绝对路径，不再假设相对路径

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(reshape2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
HOME_DIR <- path.expand("~")
DOWNLOADS_DIR <- file.path(HOME_DIR, "Downloads")
setwd(PROJECT_DIR)

message("========================================")
message("V99: F4/R3/F3 绝对路径修复")
message("========================================")

# --- 工具函数 ---
safe_read_abs <- function(path, type = "csv") {
  if (!file.exists(path)) return(NULL)
  if (type == "csv") return(read.csv(path, stringsAsFactors = FALSE))
  if (type == "rds") return(readRDS(path))
  NULL
}

REPORTS <- list()

# ========================================
# MODULE 1: F4 NYHA剂量反应（绝对路径）
# ========================================
message("")
message(">>> [F4] NYHA剂量反应修复")

tryCatch({
  # 扫描所有可能位置（Project根目录 + Downloads）
  search_paths <- c(PROJECT_DIR, DOWNLOADS_DIR)
  
  # 找V71/V78文件
  v71_file <- NULL
  v78_file <- NULL
  expr_file <- NULL
  
  for (d in search_paths) {
    if (is.null(v71_file)) {
      v71_file <- list.files(d, pattern = "V71_GSE135055_full_phenotype.*\\.csv$", full.names = TRUE)[1]
      if (is.na(v71_file)) v71_file <- NULL
    }
    if (is.null(v78_file)) {
      v78_file <- list.files(d, pattern = "V78_GSE135055.*\\.csv$", full.names = TRUE)[1]
      if (is.na(v78_file)) v78_file <- NULL
    }
    if (is.null(expr_file)) {
      expr_file <- list.files(d, pattern = "GSE135055_FPKM_Expression_Matrix.*\\.txt\\.gz$", full.names = TRUE)[1]
      if (is.na(expr_file)) expr_file <- NULL
    }
  }
  
  message("  V71表型: ", ifelse(is.null(v71_file), "未找到", basename(v71_file)))
  message("  V78临床: ", ifelse(is.null(v78_file), "未找到", basename(v78_file)))
  message("  表达矩阵: ", ifelse(is.null(expr_file), "未找到", basename(expr_file)))
  
  # 策略A: 如果V78同时有NYHA和NDUFB7，直接用
  if (!is.null(v78_file)) {
    v78 <- safe_read_abs(v78_file)
    message("  V78列名: ", paste(colnames(v78), collapse = ", "))
    
    has_nyha <- any(grepl("NYHA|nyha|class", colnames(v78), ignore.case = TRUE))
    has_ndufb7 <- any(grepl("NDUFB7|ndufb7", colnames(v78), ignore.case = TRUE))
    
    if (has_nyha && has_ndufb7) {
      message("  [FOUND] V78含NYHA+NDUFB7，直接使用")
      
      nyha_col <- grep("NYHA|nyha|class", colnames(v78), ignore.case = TRUE, value = TRUE)[1]
      ndufb7_col <- grep("NDUFB7|ndufb7", colnames(v78), ignore.case = TRUE, value = TRUE)[1]
      
      nyha_vals <- v78[[nyha_col]]
      ndufb7_vals <- v78[[ndufb7_col]]
      
      # 清理
      valid <- !is.na(nyha_vals) & !is.na(ndufb7_vals) & nyha_vals != ""
      nyha_clean <- suppressWarnings(as.numeric(gsub("[^0-9]", "", as.character(nyha_vals[valid]))))
      ndufb7_clean <- as.numeric(ndufb7_vals[valid])
      
      valid2 <- !is.na(nyha_clean) & nyha_clean >= 1 & nyha_clean <= 4 & !is.na(ndufb7_clean)
      
      if (sum(valid2) >= 5) {
        nyha_final <- nyha_clean[valid2]
        ndufb7_final <- ndufb7_clean[valid2]
        
        sp <- cor.test(ndufb7_final, nyha_final, method = "spearman")
        
        verdict <- ifelse(sp$estimate < -0.2 && sp$p.value < 0.05, "STRONG_SUPPORT",
                   ifelse(sp$p.value < 0.05, "MODERATE_SUPPORT", "NO_SUPPORT"))
        
        message("\n  [", verdict, "] NYHA vs NDUFB7:")
        message("    Spearman ρ = ", round(sp$estimate, 3), " p = ", format(sp$p.value, digits = 2))
        message("    N = ", sum(valid2))
        
        outdir <- file.path(PROJECT_DIR, "03_results/V99_F4_NYHA")
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        
        res_df <- data.frame(NYHA = nyha_final, NDUFB7 = ndufb7_final)
        write.csv(res_df, file.path(outdir, "V99_NYHA_NDUFB7.csv"), row.names = FALSE)
        
        p <- ggplot(res_df, aes(x = factor(NYHA), y = NDUFB7, fill = factor(NYHA))) +
          geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
          labs(title = "NDUFB7 by NYHA Class", 
               subtitle = paste0("ρ=", round(sp$estimate, 3), ", p=", format(sp$p.value, digits = 1))) +
          theme_minimal(base_size = 10) + theme(legend.position = "none")
        ggsave(file.path(outdir, "V99_NYHA_boxplot.png"), p, width = 6, height = 4, dpi = 300)
        
        REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = verdict, 
                                             Evidence = paste0("ρ=", round(sp$estimate, 2)), N = sum(valid2))
        message("  [DONE] F4完成 -> ", outdir)
      } else {
        message("  [FAIL] V78有效样本不足")
      }
    }
  }
  
  # 策略B: 合并V71表型 + 表达矩阵
  if (is.null(REPORTS[[length(REPORTS)]]) || REPORTS[[length(REPORTS)]]$Module != "F4") {
    if (!is.null(v71_file) && !is.null(expr_file)) {
      message("  [INFO] 尝试合并V71表型 + 表达矩阵...")
      
      pheno <- safe_read_abs(v71_file)
      expr <- safe_read_abs(expr_file)
      
      message("    V71列: ", paste(head(colnames(pheno), 8), collapse = ", "))
      message("    Expr维度: ", nrow(expr), " × ", ncol(expr))
      
      # 找基因名列（通常是第一列）
      gene_col <- colnames(expr)[1]
      ndufb7_row <- grep("NDUFB7", expr[[gene_col]], ignore.case = TRUE)
      
      if (length(ndufb7_row) > 0) {
        ndufb7_vals <- as.numeric(expr[ndufb7_row[1], -1])
        sample_ids <- colnames(expr)[-1]
        
        # 在pheno中找样本ID和NYHA
        id_col <- grep("sample|gsm|title|accession", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
        nyha_col <- grep("NYHA|nyha|class|stage|functional", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
        
        message("    匹配列: ID='", id_col, "' NYHA='", nyha_col, "'")
        
        if (!is.na(id_col) && !is.na(nyha_col)) {
          pheno_ids <- as.character(pheno[[id_col]])
          
          # 尝试匹配
          matched_indices <- sapply(sample_ids, function(s) {
            idx <- grep(s, pheno_ids, ignore.case = TRUE)[1]
            ifelse(is.na(idx), NA, idx)
          })
          
          valid <- !is.na(matched_indices)
          if (sum(valid) >= 5) {
            nyha_raw <- pheno[[nyha_col]][matched_indices[valid]]
            ndufb7_match <- ndufb7_vals[valid]
            
            # 清理NYHA
            nyha_num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", as.character(nyha_raw))))
            valid2 <- !is.na(nyha_num) & nyha_num >= 1 & nyha_num <= 4 & !is.na(ndufb7_match)
            
            if (sum(valid2) >= 5) {
              nyha_final <- nyha_num[valid2]
              ndufb7_final <- ndufb7_match[valid2]
              
              sp <- cor.test(ndufb7_final, nyha_final, method = "spearman")
              verdict <- ifelse(sp$estimate < -0.2 && sp$p.value < 0.05, "STRONG_SUPPORT",
                         ifelse(sp$p.value < 0.05, "MODERATE_SUPPORT", "NO_SUPPORT"))
              
              message("\n  [", verdict, "] 合并分析 ρ=", round(sp$estimate, 3), " p=", format(sp$p.value, digits = 2))
              
              outdir <- file.path(PROJECT_DIR, "03_results/V99_F4_NYHA")
              dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
              write.csv(data.frame(NYHA = nyha_final, NDUFB7 = ndufb7_final),
                        file.path(outdir, "V99_NYHA_NDUFB7_merged.csv"), row.names = FALSE)
              
              REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = verdict,
                                                     Evidence = "merged V71+expr", N = sum(valid2))
              message("  [DONE] F4合并完成")
            } else {
              message("  [FAIL] 合并后NYHA解析失败")
            }
          } else {
            message("  [FAIL] 样本匹配不足 (", sum(valid), ")")
          }
        } else {
          message("  [FAIL] 找不到ID/NYHA列")
        }
      } else {
        message("  [FAIL] NDUFB7不在表达矩阵")
      }
    } else {
      message("  [FAIL] V71或表达矩阵缺失")
    }
  }
  
  if (length(REPORTS) == 0 || REPORTS[[length(REPORTS)]]$Module != "F4") {
    REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "BLOCKED", Evidence = "所有策略失败", N = 0)
  }
  
}, error = function(e) {
  message("  [ERROR] F4: ", conditionMessage(e))
  REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "ERROR", Evidence = conditionMessage(e), N = 0)
})

# ========================================
# MODULE 2: R3 LVAD恢复（绝对路径+内存保护）
# ========================================
message("")
message(">>> [R3] LVAD恢复验证（绝对路径）")

tryCatch({
  # 扫描所有可能位置
  candidates <- c(
    file.path(PROJECT_DIR, "01_data/01_raw_geo/GSE226314_global.rds.gz"),
    file.path(DOWNLOADS_DIR, "GSE226314_global.rds.gz"),
    file.path(PROJECT_DIR, "Downloads/GSE226314_global.rds.gz"),
    file.path(HOME_DIR, "GSE226314_global.rds.gz")
  )
  
  global_file <- NULL
  for (c in candidates) {
    if (file.exists(c)) {
      global_file <- c
      break
    }
  }
  
  if (is.null(global_file)) {
    message("  [MISS] GSE226314_global.rds.gz 在所有位置均未找到")
    message("  搜索位置:")
    for (c in candidates) message("    - ", c)
    REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "文件缺失", N = 0)
  } else {
    sz <- file.info(global_file)$size
    message("  [FOUND] ", basename(global_file), " (", round(sz/1e9, 1), " GB)")
    message("  [INFO] 完整路径: ", global_file)
    
    # 内存检查
    mem_gb <- as.numeric(system("free -g | awk '/^Mem:/{print $2}'", intern = TRUE))
    message("  [INFO] 系统内存: ", mem_gb, " GB")
    
    if (mem_gb < 16 && sz > 5e9) {
      message("  [WARN] 内存可能不足（", mem_gb, "GB < 建议16GB），尝试加载...")
      message("  [TIP] 如果失败，建议:")
      message("    1. 关闭其他内存占用程序")
      message("    2. 或在服务器上执行")
    }
    
    message("  [LOAD] 读取中，请等待（6.4GB文件，可能需要2-5分钟）...")
    
    con <- gzfile(global_file, "rb")
    obj <- readRDS(con)
    close(con)
    
    message("  [PASS] 加载成功，类型: ", paste(class(obj), collapse = ", "))
    
    if ("Seurat" %in% class(obj)) {
      meta <- obj@meta.data
      message("  Meta列: ", paste(head(colnames(meta), 10), collapse = ", "))
      
      # 找分组列
      grp_candidates <- c("condition", "Condition", "time_point", "treatment", "LVAD", "group", "status", "disease_state")
      grp_col <- intersect(grp_candidates, colnames(meta))[1]
      
      if (!is.na(grp_col)) {
        grps <- unique(meta[[grp_col]])
        message("  分组列 '", grp_col, "': ", paste(grps, collapse = ", "))
        
        pre_mask <- grepl("pre|baseline|before|0h|0$|prior", meta[[grp_col]], ignore.case = TRUE)
        post_mask <- grepl("post|after|LVAD|1w|4w|follow", meta[[grp_col]], ignore.case = TRUE)
        
        message("  Pre匹配: ", sum(pre_mask), " cells")
        message("  Post匹配: ", sum(post_mask), " cells")
        
        if (sum(pre_mask) > 5 && sum(post_mask) > 5) {
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
            
            outdir <- file.path(PROJECT_DIR, "03_results/V99_R3_LVAD")
            dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
            
            res_df <- data.frame(
              Group = c(rep("Pre-LVAD", length(pre_vals)), rep("Post-LVAD", length(post_vals))),
              NDUFB7 = c(pre_vals, post_vals)
            )
            write.csv(res_df, file.path(outdir, "V99_LVAD_NDUFB7.csv"), row.names = FALSE)
            
            p <- ggplot(res_df, aes(x = Group, y = NDUFB7, fill = Group)) +
              geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
              annotate("text", x = 1.5, y = Inf,
                       label = paste0("Δ=", round(pct, 1), "%, p=", format(tt$p.value, digits = 1)),
                       vjust = 1.5, size = 3) +
              labs(title = "NDUFB7 Recovery After LVAD") +
              theme_minimal(base_size = 10) + theme(legend.position = "none")
            ggsave(file.path(outdir, "V99_LVAD_boxplot.png"), p, width = 5, height = 4, dpi = 300)
            
            REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = verdict,
                                                 Evidence = paste0("Δ=", round(pct, 1), "%"), 
                                                 N = paste0(length(pre_vals), "/", length(post_vals)))
            message("  [DONE] R3完成 -> ", outdir)
            
          } else {
            message("  [FAIL] NDUFB7不在基因列表")
            REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "NDUFB7缺失", N = 0)
          }
        } else {
          message("  [FAIL] Pre/Post分组不足")
          REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "分组不足", N = 0)
        }
      } else {
        message("  [FAIL] 无分组列")
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
    message("  [DIAGNOSIS] 内存不足")
    message("  [WORKAROUND] 在Discussion写: 'Future validation in independent LVAD cohort required'")
  }
  REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "ERROR", Evidence = msg, N = 0)
})

# ========================================
# MODULE 3: F3 GSE243655 Fer-1（绝对路径归档+解析）
# ========================================
message("")
message(">>> [F3] GSE243655 Fer-1数据准备")

tryCatch({
  # 扫描tar文件
  tar_candidates <- c(
    file.path(PROJECT_DIR, "01_data/01_raw_geo/GSE243655/GSE243655_RAW(1).tar"),
    file.path(DOWNLOADS_DIR, "GSE243655_RAW(1).tar"),
    file.path(DOWNLOADS_DIR, "GSE243655_RAW.tar"),
    file.path(PROJECT_DIR, "Downloads/GSE243655_RAW(1).tar")
  )
  
  tar_file <- NULL
  for (c in tar_candidates) {
    if (file.exists(c)) {
      tar_file <- c
      break
    }
  }
  
  if (is.null(tar_file)) {
    message("  [MISS] GSE243655 RAW.tar 未找到")
    message("  搜索位置:")
    for (c in tar_candidates) message("    - ", c)
    REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "BLOCKED", Evidence = "文件缺失", N = 0)
  } else {
    message("  [FOUND] ", basename(tar_file))
    message("  完整路径: ", tar_file)
    
    # 解压到Project标准目录
    EXTRACT_DIR <- file.path(PROJECT_DIR, "01_data/01_raw_geo/GSE243655")
    dir.create(EXTRACT_DIR, showWarnings = FALSE, recursive = TRUE)
    
    existing <- list.files(EXTRACT_DIR, pattern = "\\.txt$|\\.csv$|\\.cel$", ignore.case = TRUE)
    
    if (length(existing) == 0) {
      message("  [EXTRACT] 解压到 ", EXTRACT_DIR)
      system(paste("tar -xf", shQuote(tar_file), "-C", shQuote(EXTRACT_DIR)), intern = TRUE)
      existing <- list.files(EXTRACT_DIR, pattern = "\\.txt$|\\.csv$|\\.cel$", ignore.case = TRUE)
    }
    
    message("  [INFO] 文件数: ", length(existing))
    message("  文件列表:")
    for (f in head(existing, 10)) message("    - ", f)
    
    # 尝试读取series matrix
    sm_candidates <- c(
      file.path(EXTRACT_DIR, "GSE243655_series_matrix.txt.gz"),
      file.path(DOWNLOADS_DIR, "GSE243655_series_matrix.txt.gz"),
      file.path(PROJECT_DIR, "GSE243655_series_matrix.txt.gz")
    )
    sm_file <- NULL
    for (c in sm_candidates) {
      if (file.exists(c)) {
        sm_file <- c
        break
      }
    }
    
    if (!is.null(sm_file)) {
      message("  [INFO] Series matrix: ", basename(sm_file))
      # 读取并提取分组
      sm_lines <- readLines(sm_file, n = 200)
      char_lines <- sm_lines[grep("!Sample_characteristics", sm_lines)]
      
      message("  特征行:")
      for (l in head(char_lines, 5)) message("    ", substr(l, 1, 80))
      
      # 提取treatment信息
      treat_lines <- sm_lines[grep("treatment|agent|drug|ferrostatin|vehicle", sm_lines, ignore.case = TRUE)]
      message("  Treatment行:")
      for (l in head(treat_lines, 5)) message("    ", substr(l, 1, 80))
      
      # 保存分组信息供手动解析
      outdir <- file.path(PROJECT_DIR, "03_results/V99_F3_Fer1")
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
      writeLines(sm_lines, file.path(outdir, "V99_GSE243655_series_matrix_header.txt"))
      
      message("\n  [PENDING_MANUAL] 数据已解压，分组信息已保存")
      message("  [ACTION] 请查看: ", file.path(outdir, "V99_GSE243655_series_matrix_header.txt"))
      message("  [ACTION] 确认Fer-1 vs Vehicle分组后，执行后续差异分析")
      
      REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "PENDING_MANUAL",
                                           Evidence = paste0("解压完成，", length(existing), "个文件"),
                                           N = length(existing))
    } else {
      message("  [INFO] 无series matrix，需手动确认分组")
      REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "PENDING_MANUAL",
                                           Evidence = "无series matrix", N = length(existing))
    }
    message("  [DONE] F3数据准备完成")
  }
  
}, error = function(e) {
  message("  [ERROR] F3: ", conditionMessage(e))
  REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "ERROR", Evidence = conditionMessage(e), N = 0)
})

# ========================================
# 汇总
# ========================================
message("")
message("========================================")
message("V99 修复结果汇总")
message("========================================")

report_df <- do.call(rbind, lapply(REPORTS, function(x) {
  data.frame(Module = x$Module, Verdict = x$Verdict, Evidence = x$Evidence, N = as.character(x$N), stringsAsFactors = FALSE)
}))
print(report_df, row.names = FALSE)

outdir <- file.path(PROJECT_DIR, "03_results/V99_Validation_Fix")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(report_df, file.path(outdir, "V99_fix_report.csv"), row.names = FALSE)

message("")
message("[DONE] V99完成")
message("  R4(已修复): 03_results/V98_R4_Death_Independence/ (铁死亡独立性✅)")
message("  F4: 03_results/V99_F4_NYHA/")
message("  R3: 03_results/V99_R3_LVAD/")
message("  F3: 01_data/01_raw_geo/GSE243655/ + 03_results/V99_F3_Fer1/")
