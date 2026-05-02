#!/usr/bin/env Rscript
# V100: 终极修复——硬编码V97已知路径 + 多格式探测
# 不再依赖list.files，直接指向V97诊断确认的文件

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

PROJECT_DIR <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
HOME_DIR <- path.expand("~")
setwd(PROJECT_DIR)

message("========================================")
message("V100: 终极修复（硬编码路径 + 格式探测）")
message("========================================")

REPORTS <- list()

# ========================================
# MODULE 1: F4 NYHA（硬编码V97确认的路径）
# ========================================
message("")
message(">>> [F4] NYHA剂量反应（硬编码路径）")

tryCatch({
  # V97诊断确认的文件路径（直接硬编码，不搜索）
  V78_PATH <- file.path(PROJECT_DIR, "V78_GSE135055_NDUFB7_clinical.csv")
  V71_PATH <- file.path(PROJECT_DIR, "V71_GSE135055_full_phenotype_SOFT.csv")
  V71_RNA_PATH <- file.path(PROJECT_DIR, "V71_GSE135055_RNAseq_phenotype.csv")
  EXPR_PATH <- file.path(HOME_DIR, "Downloads", "GSE135055_FPKM_Expression_Matrix(1).txt.gz")
  EXPR_PATH2 <- file.path(HOME_DIR, "Downloads", "GSE135055_FPKM_Expression_Matrix.txt.gz")
  
  # 策略A: 直接读取V78（如果含NYHA+NDUFB7）
  if (file.exists(V78_PATH)) {
    message("  [FOUND] V78: ", basename(V78_PATH))
    v78 <- read.csv(V78_PATH, stringsAsFactors = FALSE)
    message("  列名: ", paste(colnames(v78), collapse = ", "))
    message("  行数: ", nrow(v78))
    
    # 检查列
    has_nyha <- any(grepl("NYHA|nyha|class|stage", colnames(v78), ignore.case = TRUE))
    has_ndufb7 <- any(grepl("NDUFB7|ndufb7", colnames(v78), ignore.case = TRUE))
    
    if (has_nyha && has_ndufb7) {
      message("  [PASS] V78含NYHA+NDUFB7，直接分析")
      
      nyha_col <- grep("NYHA|nyha|class|stage", colnames(v78), ignore.case = TRUE, value = TRUE)[1]
      ndufb7_col <- grep("NDUFB7|ndufb7", colnames(v78), ignore.case = TRUE, value = TRUE)[1]
      
      nyha_raw <- v78[[nyha_col]]
      ndufb7_raw <- suppressWarnings(as.numeric(v78[[ndufb7_col]]))
      
      # 清理NYHA
      nyha_num <- suppressWarnings(as.numeric(gsub("[^0-9]", "", as.character(nyha_raw))))
      valid <- !is.na(nyha_num) & nyha_num >= 1 & nyha_num <= 4 & !is.na(ndufb7_raw)
      
      if (sum(valid) >= 5) {
        nyha_final <- nyha_num[valid]
        ndufb7_final <- ndufb7_raw[valid]
        
        sp <- cor.test(ndufb7_final, nyha_final, method = "spearman")
        
        verdict <- ifelse(sp$estimate < -0.2 && sp$p.value < 0.05, "STRONG_SUPPORT",
                   ifelse(sp$p.value < 0.05, "MODERATE_SUPPORT", "NO_SUPPORT"))
        
        message("\n  [", verdict, "]")
        message("    Spearman ρ = ", round(sp$estimate, 3), " p = ", format(sp$p.value, digits = 2))
        message("    N = ", sum(valid))
        
        outdir <- file.path(PROJECT_DIR, "03_results/V100_F4_NYHA")
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        
        res_df <- data.frame(NYHA = nyha_final, NDUFB7 = ndufb7_final)
        write.csv(res_df, file.path(outdir, "V100_NYHA_NDUFB7.csv"), row.names = FALSE)
        
        p <- ggplot(res_df, aes(x = factor(NYHA), y = NDUFB7, fill = factor(NYHA))) +
          geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
          labs(title = "NDUFB7 by NYHA Class",
               subtitle = paste0("ρ=", round(sp$estimate, 3), ", p=", format(sp$p.value, digits = 1))) +
          theme_minimal(base_size = 10) + theme(legend.position = "none")
        ggsave(file.path(outdir, "V100_NYHA_boxplot.png"), p, width = 6, height = 4, dpi = 300)
        
        REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = verdict,
                                             Evidence = paste0("V78 ρ=", round(sp$estimate, 2)), N = sum(valid))
        message("  [DONE] F4完成 -> ", outdir)
      } else {
        message("  [FAIL] V78有效样本不足 (", sum(valid), ")")
      }
    } else {
      message("  [INFO] V78缺少NYHA或NDUFB7，尝试V71...")
    }
  }
  
  # 策略B: V71 + 表达矩阵合并
  if (length(REPORTS) == 0 || REPORTS[[length(REPORTS)]]$Module != "F4") {
    if (file.exists(V71_PATH) && (file.exists(EXPR_PATH) || file.exists(EXPR_PATH2))) {
      expr_path <- ifelse(file.exists(EXPR_PATH), EXPR_PATH, EXPR_PATH2)
      message("  [FOUND] V71 + 表达矩阵")
      
      pheno <- read.csv(V71_PATH, stringsAsFactors = FALSE)
      message("    V71列: ", paste(head(colnames(pheno), 8), collapse = ", "))
      
      # 尝试读取表达矩阵（大文件，可能慢）
      message("    [LOAD] 读取表达矩阵...")
      expr <- read.table(expr_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, nrows = 5)
      message("    Expr预览: ", nrow(expr), " × ", ncol(expr))
      
      # 这里简化：如果V71本身有NDUFB7和NYHA，直接用
      # 否则标记为需要手动合并
      message("  [PENDING] V71+表达矩阵合并需要手动匹配样本ID")
      REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "PENDING_MANUAL",
                                           Evidence = "V71和表达矩阵存在，需手动匹配样本ID", N = 0)
    } else {
      message("  [FAIL] V71或表达矩阵不存在")
      message("    V71存在: ", file.exists(V71_PATH))
      message("    Expr存在: ", file.exists(EXPR_PATH), " / ", file.exists(EXPR_PATH2))
      REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "BLOCKED", Evidence = "文件缺失", N = 0)
    }
  }
  
}, error = function(e) {
  message("  [ERROR] F4: ", conditionMessage(e))
  REPORTS[[length(REPORTS)+1]] <- list(Module = "F4", Verdict = "ERROR", Evidence = conditionMessage(e), N = 0)
})

# ========================================
# MODULE 2: R3 LVAD（格式探测 + 多方式读取）
# ========================================
message("")
message(">>> [R3] LVAD恢复验证（格式探测）")

tryCatch({
  R3_PATH <- file.path(HOME_DIR, "Downloads", "GSE226314_global.rds.gz")
  
  if (!file.exists(R3_PATH)) {
    message("  [MISS] 文件不存在")
    REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "文件缺失", N = 0)
  } else {
    sz <- file.info(R3_PATH)$size
    message("  [FOUND] ", basename(R3_PATH), " (", round(sz/1e9, 1), " GB)")
    
    # 格式探测：读取前100字节
    con <- gzfile(R3_PATH, "rb")
    header <- readBin(con, "raw", n = 100)
    close(con)
    header_str <- rawToChar(header[header != 0])
    
    message("  [PROBE] 文件头: ", substr(header_str, 1, 50))
    
    # 判断格式
    is_hdf5 <- any(header[1:8] == as.raw(c(0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a)))
    is_rds <- grepl("RDS|RDA|RDATA|A\\x02|B\\x02", header_str)
    
    message("  HDF5签名: ", is_hdf5)
    message("  RDS签名: ", is_rds)
    
    if (is_hdf5) {
      message("  [FORMAT] 检测到HDF5格式（可能是h5Seurat或h5ad）")
      message("  [FIX] 尝试用SeuratDisk读取...")
      
      if (requireNamespace("SeuratDisk", quietly = TRUE)) {
        library(SeuratDisk)
        obj <- LoadH5Seurat(R3_PATH)
        message("  [PASS] h5Seurat加载成功")
        # 后续分析同前...
      } else {
        message("  [FAIL] SeuratDisk未安装")
        message("  [ACTION] install.packages('SeuratDisk') 或 conda install -c conda-forge r-seuratdisk")
        REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "HDF5格式需SeuratDisk", N = 0)
      }
    } else if (is_rds) {
      message("  [FORMAT] 检测到RDS格式，但可能损坏或版本不兼容")
      message("  [ACTION] 尝试在R 4.3+中读取，或检查文件完整性")
      REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "RDS格式损坏或版本不兼容", N = 0)
    } else {
      message("  [FORMAT] 未知格式，可能为自定义压缩或加密")
      message("  [ACTION] 联系数据提供者确认格式，或尝试用Python scanpy读取")
      REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "BLOCKED", Evidence = "未知文件格式", N = 0)
    }
  }
  
}, error = function(e) {
  message("  [ERROR] R3: ", conditionMessage(e))
  REPORTS[[length(REPORTS)+1]] <- list(Module = "R3", Verdict = "ERROR", Evidence = conditionMessage(e), N = 0)
})

# ========================================
# MODULE 3: F3 Fer-1（解析sample_list + 指导GEO下载）
# ========================================
message("")
message(">>> [F3] GSE243655 Fer-1分组解析")

tryCatch({
  EXTRACT_DIR <- file.path(PROJECT_DIR, "01_data/01_raw_geo/GSE243655")
  SAMPLE_LIST <- file.path(EXTRACT_DIR, "sample_list.txt")
  SM_PATH <- file.path(HOME_DIR, "Downloads", "GSE243655_series_matrix.txt.gz")
  
  # 读取sample_list
  if (file.exists(SAMPLE_LIST)) {
    samples <- readLines(SAMPLE_LIST)
    message("  Sample list (", length(samples), " lines):")
    for (s in head(samples, 10)) message("    ", s)
  }
  
  # 解析series matrix提取分组
  if (file.exists(SM_PATH)) {
    message("  [FOUND] Series matrix")
    sm_lines <- readLines(SM_PATH, n = 300)
    
    # 提取所有!Sample_characteristics_ch1行
    char_lines <- sm_lines[grep("^!Sample_characteristics_ch1", sm_lines)]
    treat_lines <- sm_lines[grep("treatment|DMSO|ferrostatin|vehicle", sm_lines, ignore.case = TRUE)]
    
    message("\n  === 关键特征行 ===")
    for (l in head(char_lines, 6)) {
      clean <- gsub("!Sample_characteristics_ch1\t", "", l)
      message("    ", substr(clean, 1, 100))
    }
    
    # 提取分组
    # 寻找"treatment:"开头的行
    treat_idx <- grep("treatment:", char_lines, ignore.case = TRUE)
    if (length(treat_idx) > 0) {
      treat_line <- char_lines[treat_idx[1]]
      # 解析tab分隔的值
      treat_vals <- strsplit(treat_line, "\t")[[1]][-1]  # 去掉标题
      treat_vals <- gsub("^\"|\"$", "", treat_vals)     # 去掉引号
      
      message("\n  === Treatment分组 ===")
      tbl <- table(treat_vals)
      for (nm in names(tbl)) {
        message("    ", nm, ": ", tbl[nm], " samples")
      }
      
      # 判断是否有Fer-1 vs DMSO/Vehicle
      has_ferro <- any(grepl("ferrostatin", treat_vals, ignore.case = TRUE))
      has_dmso <- any(grepl("DMSO|vehicle|control", treat_vals, ignore.case = TRUE))
      
      if (has_ferro && has_dmso) {
        message("\n  [PASS] 确认分组: Fer-1 vs DMSO/Vehicle")
        message("  [ACTION] 需要从GEO下载表达矩阵，然后:")
        message("    1. 提取NDUFB7表达值")
        message("    2. 按treatment分组做t-test")
        message("    3. 比较Fer-1 vs DMSO的NDUFB7变化")
        
        # 保存分组表
        outdir <- file.path(PROJECT_DIR, "03_results/V100_F3_Fer1")
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        
        grp_df <- data.frame(
          Sample = seq_along(treat_vals),
          Treatment = treat_vals,
          Group = ifelse(grepl("ferrostatin", treat_vals, ignore.case = TRUE), "Ferrostatin-1",
                  ifelse(grepl("DMSO|vehicle", treat_vals, ignore.case = TRUE), "Vehicle", "Other"))
        )
        write.csv(grp_df, file.path(outdir, "V100_GSE243655_groups.csv"), row.names = FALSE)
        
        REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "PENDING_EXPRESSION",
                                             Evidence = paste0("Fer-1=", sum(grepl("ferrostatin", treat_vals, ignore.case = TRUE)),
                                                              " DMSO=", sum(grepl("DMSO", treat_vals, ignore.case = TRUE))),
                                             N = length(treat_vals))
        message("  [DONE] 分组表保存 -> ", outdir)
      } else {
        message("  [WARN] 未识别Fer-1 vs DMSO分组")
        REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "PENDING_MANUAL", Evidence = "分组不明确", N = 0)
      }
    } else {
      message("  [FAIL] Series matrix中无treatment行")
      REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "BLOCKED", Evidence = "无treatment信息", N = 0)
    }
  } else {
    message("  [MISS] Series matrix不存在")
    REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "BLOCKED", Evidence = "无series matrix", N = 0)
  }
  
}, error = function(e) {
  message("  [ERROR] F3: ", conditionMessage(e))
  REPORTS[[length(REPORTS)+1]] <- list(Module = "F3", Verdict = "ERROR", Evidence = conditionMessage(e), N = 0)
})

# ========================================
# 汇总 + 论文写作建议
# ========================================
message("")
message("========================================")
message("V100 修复结果汇总")
message("========================================")

report_df <- do.call(rbind, lapply(REPORTS, function(x) {
  data.frame(Module = x$Module, Verdict = x$Verdict, Evidence = x$Evidence, N = as.character(x$N), stringsAsFactors = FALSE)
}))
print(report_df, row.names = FALSE)

outdir <- file.path(PROJECT_DIR, "03_results/V100_Validation_Final")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(report_df, file.path(outdir, "V100_fix_report.csv"), row.names = FALSE)

# 生成论文建议
message("")
message("========================================")
message("假说评估与论文写作建议")
message("========================================")

# 统计
pass_fwd <- sum(report_df$Verdict %in% c("STRONG_SUPPORT", "MODERATE_SUPPORT", "PASS", "INDEPENDENT"))
pass_rev <- sum(report_df$Verdict %in% c("REVERSED", "PARTIALLY_REVERSED", "INDEPENDENT"))
pending <- sum(report_df$Verdict %in% c("PENDING_MANUAL", "PENDING_EXPRESSION", "BLOCKED", "ERROR"))

message("正向验证通过: ", pass_fwd, " 项")
message("反向验证通过: ", pass_rev, " 项")
message("待处理/阻塞:  ", pending, " 项")

message("")
message("## 核心结论")

if (pass_fwd >= 2 && pass_rev >= 1 && pending <= 2) {
  message("[GRADE A] 假说可通过压力测试，建议写入论文")
  message("  核心叙事: NDUFB7↓作为OXPHOS崩溃-铁死亡轴的关键节点")
  message("  可claim: '关联性' + '机制推断'，避免直接claim 'causation'")
} else if (pass_fwd >= 2 && pass_rev >= 1) {
  message("[GRADE B] 假说部分支持，需弱化表述")
  message("  核心叙事: NDUFB7↓与OXPHOS崩溃、铁死亡激活形成病理三联征")
  message("  必须避免: 'NDUFB7驱动HF'、'铁死亡是核心机制'等强因果表述")
  message("  建议增加: '本研究为观察性关联，功能实验验证有待后续研究'")
} else {
  message("[GRADE C] 假说证据不足，建议大幅降级")
  message("  核心叙事: NDUFB7是HF的潜在生物标志物")
  message("  删除所有机制推断，仅保留差异表达和预后关联")
}

message("")
message("## 具体写作策略")

message("\n### Results 可写段落（有数据支撑）:")
message("1. Fig 1: NDUFB7在HF中下调（GSE183852 + meta-analysis）")
message("2. Fig 2: CM亚群分层 + HOX亚群鉴定（V81 + V93）")
message("3. Fig 3: OXPHOS多复合体崩溃（T27, I/III/IV下调）")
message("4. Fig 4: 铁死亡独立性验证（R4, ρ≈0）")
message("5. Fig 5: 空间异质性（Visium, 如Moran's I重建成功）")
message("6. Fig 6: MR因果推断（已有eQTLGen/GTEx/HERMES）")

message("\n### Discussion 可写段落（机制推测，带保守措辞）:")
message("- 'Our findings suggest that NDUFB7 loss may contribute to a vicious cycle")
message("  of mitochondrial dysfunction and ferroptosis...'（用may, suggest）")
message("- 'Future studies using genetic knockout or Fer-1 intervention are needed")
message("  to establish causal relationships...'（主动承认局限）")

message("\n### 必须删除/避免:")
message("- 'NDUFB7 drives HF through ferroptosis'（无Fer-1挽救证据）")
message("- 'Therapeutic targeting of NDUFB7'（无LVAD恢复或药物验证）")
message("- 'NDUFB7 is the core mechanism'（无细胞特异性证据）")

message("\n### 审稿人可能质疑及回应:")
message("Q: 'How do you rule out apoptosis as the real cause of CM death?'")
message("A: 'We showed ferroptosis, apoptosis, and pyroptosis are nearly")
message("   uncorrelated (R4, ρ=0.02), suggesting independent pathways.'")

message("Q: 'Is NDUFB7 loss reversible? If not, why target it?'")
message("A: 'While our LVAD analysis was inconclusive due to data format issues,")
message("   the compensatory upregulation in HOX cells (T3) suggests dynamic")
message("   regulation, supporting therapeutic relevance.'")

message("[DONE] 完整报告: ", outdir)
