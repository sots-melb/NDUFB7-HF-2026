#!/usr/bin/env Rscript
PROJECT_ROOT <- "~/Projects/NDUFB7_HF_{2026_04_20}"
setwd(PROJECT_ROOT)

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║     GEO 2万+摘要解析 & 空间转录组智能筛选                ║\n")
cat("║     设计原则: 高效、可审计、分级输出、一键下游           ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")

# ============================================================
# 0. 自动探测GEO摘要文件
# ============================================================
cat("\n========== 0. 探测GEO Summary文件 ==========\n")

search_dirs <- c(
  "01_data/01_raw_geo/",
  "~/Downloads/",
  "~/download/",
  "~/Desktop/",
  ".",
  "../"
)

possible_patterns <- c(
  "geo.*summary.*\\.txt",
  "summary.*geo.*\\.txt",
  "heart.*mi.*summary",
  "ncbi.*geo.*summary",
  "geo.*series.*summary",
  "gse.*summary",
  ".*summary.*\\.txt$"
)

summary_file <- NULL
for (dir in search_dirs) {
  dir_exp <- path.expand(dir)
  if (!dir.exists(dir_exp)) next
  
  files <- list.files(dir_exp, pattern = "\\.txt$|\\.csv$", full.names = TRUE, ignore.case = TRUE)
  for (f in files) {
    # 快速检查：文件大小>10KB且包含GSE
    f_info <- file.info(f)
    if (f_info$size > 10000) {
      first_lines <- readLines(f, n = 50, warn = FALSE)
      if (any(grepl("GSE\\d+", first_lines))) {
        summary_file <- f
        cat("✅ 找到GEO摘要文件:\n")
        cat("   路径:", f, "\n")
        cat("   大小:", round(f_info$size/1024/1024, 2), "MB (", f_info$size, "bytes )\n")
        cat("   前3行预览:\n")
        cat("   ", head(first_lines, 3), sep = "\n   ")
        cat("\n")
        break
      }
    }
  }
  if (!is.null(summary_file)) break
}

if (is.null(summary_file)) {
  cat("\n❌ 未自动找到GEO摘要文件。\n")
  cat("【手动指定】请执行:\n")
  cat('   summary_file <- "您的实际路径.txt"\n')
  cat("   然后重新运行此脚本。\n")
  cat("\n【常见位置检查】:\n")
  for (dir in search_dirs) {
    dir_exp <- path.expand(dir)
    if (dir.exists(dir_exp)) {
      txt_files <- list.files(dir_exp, pattern = "\\.txt$", ignore.case = TRUE)
      if (length(txt_files) > 0) {
        cat("   ", dir, ":", paste(head(txt_files, 5), collapse = ", "), "\n")
      }
    }
  }
  quit(status = 1)
}

# ============================================================
# 1. 高效解析2万+GEO记录
# ============================================================
cat("\n========== 1. 解析GEO记录 (预计1-3分钟) ==========\n")

lines <- readLines(summary_file, warn = FALSE)
total_lines <- length(lines)
cat("总行数:", total_lines, "\n")

# --- 格式探测 ---
cat("⏳ 探测记录分隔格式...\n")

# 策略A: 编号列表 (1. GSE12345 ...)
numbered_pattern <- "^\\s*\\d+\\.[\\s\\t]+GSE\\d+"
numbered_idx <- grep(numbered_pattern, lines)

# 策略B: 纯GSE开头
gse_only_pattern <- "^\\s*GSE\\d+\\s*$"
gse_only_idx <- grep(gse_only_pattern, lines)

# 策略C: "Accession: GSE" 或 "Series: GSE"
accession_pattern <- "(?i)^(accession|series|geo accession)[:\\s]+GSE\\d+"
accession_idx <- grep(accession_pattern, lines, perl = TRUE)

cat("   编号格式候选:", length(numbered_idx), "行\n")
cat("   GSE纯名候选:", length(gse_only_idx), "行\n")
cat("   Accession候选:", length(accession_idx), "行\n")

# 选择最佳策略
if (length(numbered_idx) > 100) {
  cat("✅ 采用策略A: 编号列表格式 (", length(numbered_idx), "条记录)\n")
  record_starts <- numbered_idx
} else if (length(gse_only_idx) > 100) {
  cat("✅ 采用策略B: GSE分隔符格式 (", length(gse_only_idx), "条记录)\n")
  record_starts <- gse_only_idx
} else if (length(accession_idx) > 100) {
  cat("✅ 采用策略C: Accession格式 (", length(accession_idx), "条记录)\n")
  record_starts <- accession_idx
} else {
  # 策略D: 任何包含GSE的行，去重连续行
  all_gse <- grep("GSE\\d{3,}", lines)
  # 保留间隔>2的行（避免同一记录内多次匹配）
  record_starts <- all_gse[c(TRUE, diff(all_gse) > 2)]
  cat("✅ 采用策略D: 启发式GSE检测 (", length(record_starts), "条记录)\n")
}

if (length(record_starts) == 0) {
  cat("❌ 无法识别记录格式。请检查文件内容。\n")
  quit(status = 1)
}

cat("预计记录数:", length(record_starts), "\n")

# --- 分块提取 ---
cat("⏳ 分块提取记录文本 (向量化处理)...\n")

# 预分配
n_records <- length(record_starts)
record_texts <- character(n_records)

for (i in 1:n_records) {
  start <- record_starts[i]
  end <- ifelse(i < n_records, record_starts[i+1] - 1, total_lines)
  record_texts[i] <- paste(lines[start:end], collapse = " ")
}

# 提取GSE编号
gse_ids <- regmatches(record_texts, regexpr("GSE\\d+", record_texts))
cat("✅ 提取GSE编号完成，前5个:", paste(head(gse_ids, 5), collapse = ", "), "\n")

# --- 字段提取函数（高效正则）---
extract_field_fast <- function(texts, field_name) {
  # 匹配 "Field:" 或 "Field " 后到下一个字段或结束
  pattern <- paste0("(?i)(?:^|\\n|\\s)", field_name, "\\s*[:\\-]?\\s*(.+?)(?=\\n[A-Z][A-Za-z\\s]{2,20}[:\\-]|\\n\\d+\\.\\s|\\nGSE\\d|$)")
  result <- character(length(texts))
  
  for (i in 1:length(texts)) {
    m <- regexpr(pattern, texts[i], perl = TRUE)
    if (m[1] != -1) {
      capture <- regmatches(texts[i], m)[1]
      # 去掉字段名前缀
      clean <- sub(paste0("(?i)(?:^|\\n|\\s)", field_name, "\\s*[:\\-]?\\s*"), "", capture, perl = TRUE)
      result[i] <- trimws(clean)
    } else {
      result[i] <- NA
    }
  }
  return(result)
}

cat("⏳ 提取关键字段...\n")
fields_to_extract <- c("Title", "Summary", "Organism", "Platform", "Type", "Status", "PubDate", "Samples")
records_df <- data.frame(
  GSE = gse_ids,
  stringsAsFactors = FALSE
)

for (fld in fields_to_extract) {
  cat("   提取", fld, "...\n")
  records_df[[fld]] <- extract_field_fast(record_texts, fld)
  # 截断过长文本
  records_df[[fld]] <- substr(records_df[[fld]], 1, 500)
}

# 尝试从其他格式补全缺失的Samples
if (all(is.na(records_df$Samples))) {
  # 尝试匹配 "Samples (45)" 或 "Sample(s): 45"
  samples_pattern <- "(?i)sample[s]?\\s*\\(?s?\\)?\\s*[:\\-]?\\s*(\\d+)"
  samples_extracted <- regmatches(record_texts, regexpr(samples_pattern, record_texts, perl = TRUE))
  samples_nums <- regmatches(samples_extracted, regexpr("\\d+", samples_extracted))
  records_df$Samples <- as.character(samples_nums)
}

cat("✅ 字段提取完成。样本结构:\n")
print(head(records_df[, c("GSE", "Title", "Organism", "Type", "Samples")], 3), row.names = FALSE)

# 保存解析结果
write.csv(records_df, "03_results/16_spatial_from_bulk/01_all_parsed_records.csv", row.names = FALSE)
cat("✅ 全量解析结果已保存 (", nrow(records_df), "条)\n")

# ============================================================
# 2. 空间转录组智能筛选
# ============================================================
cat("\n========== 2. 空间转录组智能筛选 ==========\n")

# --- 关键词定义 ---
keywords <- list(
  spatial_tech = c("spatial", "visium", "slide-seq", "slide seq", "stereoseq", "stereo-seq",
                   "seq-scope", "merfish", "osmfish", "starmap", "exseq", "fisseq",
                   "10x genomics visium", "spatial transcriptomic", "spatial transcriptome",
                   "in situ sequencing", "xenium", "cosmx", "nanostring geo", "dbit-seq",
                   "hdst", "pixseq", "proseq", "seqfish"),
  
  heart = c("heart", "cardiac", "myocardial", "cardiomyocyte", "ventricle", "atrium",
            "heart failure", " hf ", "myocardial infarction", " mi ", "ischemia",
            "reperfusion", "cardiomyopathy", "dcm", "dilated cardiomyopathy",
            "hcm", "hypertrophic cardiomyopathy", "arrhythmia", "fibrosis",
            "tac", "transverse aortic constriction", "pressure overload",
            "aortic stenosis", "heart transplant", "lvad", "left ventricular"),
  
  human = c("homo sapiens", "human", "patient", "clinical", " Homo "),
  
  mouse = c("mus musculus", "mouse", "mice", "rat ", "rattus", "rodent"),
  
  scrna = c("single cell", "scrna", "single-cell", "10x chromium", "drop-seq",
            "smart-seq", "smartseq", "indrop", "seq-well"),
  
  exclude = c("methylation array", "chip-seq", "atac-seq", "cut&run", "cut and run",
              "genome sequencing", "exome sequencing", "whole genome", "mirna", "microrna",
              "methylation profiling", "dna methylation", "bisulfite", "epigenomic",
              "proteomic", "metabolomic", "genotyping array", "snp array", "cnv",
              "ribosome profiling", "rip-seq", "clip-seq", "gwas", "ewas")
)

# --- 向量化评分 ---
score_texts <- tolower(record_texts)

cat("⏳ 计算空间评分...\n")
spatial_scores <- sapply(keywords$spatial_tech, function(k) grepl(k, score_texts, fixed = TRUE))
records_df$spatial_score <- rowSums(spatial_scores) * 3  # 空间技术权重高

cat("⏳ 计算心脏评分...\n")
heart_scores <- sapply(keywords$heart, function(k) grepl(k, score_texts, fixed = TRUE))
records_df$heart_score <- rowSums(heart_scores) * 2

cat("⏳ 计算物种评分...\n")
records_df$is_human <- grepl(paste(keywords$human, collapse = "|"), score_texts, ignore.case = TRUE)
records_df$is_mouse <- grepl(paste(keywords$mouse, collapse = "|"), score_texts, ignore.case = TRUE)
records_df$species_score <- ifelse(records_df$is_human, 10, ifelse(records_df$is_mouse, 5, 2))

cat("⏳ 计算单细胞加分...\n")
scrna_scores <- sapply(keywords$scrna, function(k) grepl(k, score_texts, fixed = TRUE))
records_df$sc_bonus <- ifelse(rowSums(scrna_scores) > 0, 2, 0)

cat("⏳ 计算排除扣分...\n")
exclude_scores <- sapply(keywords$exclude, function(k) grepl(k, score_texts, fixed = TRUE))
records_df$exclude_penalty <- rowSums(exclude_scores) * 8
records_df$should_exclude <- records_df$exclude_penalty > 0

# 综合得分
records_df$total_score <- pmax(0, records_df$spatial_score + records_df$heart_score + 
                                  records_df$species_score + records_df$sc_bonus - 
                                  records_df$exclude_penalty)

# --- 分级 ---
cat("\n========== 3. 分级输出 ==========\n")

# P0: 人 + 心脏 + 空间 (最高优先级)
p0 <- records_df[records_df$spatial_score > 0 & records_df$heart_score > 0 & 
                   records_df$is_human & !records_df$should_exclude, ]
p0 <- p0[order(p0$total_score, decreasing = TRUE), ]

# P1: 人 + 心脏 + 单细胞 (空间整合备选)
p1 <- records_df[records_df$spatial_score == 0 & records_df$heart_score > 0 & 
                   records_df$is_human & records_df$sc_bonus > 0 & !records_df$should_exclude, ]
p1 <- p1[order(p1$total_score, decreasing = TRUE), ]

# P2: 小鼠 + 心脏 + 空间 (跨物种)
p2 <- records_df[records_df$spatial_score > 0 & records_df$heart_score > 0 & 
                   records_df$is_mouse & !records_df$should_exclude, ]
p2 <- p2[order(p2$total_score, decreasing = TRUE), ]

# P3: 人 + 心脏 + Bulk (扩展验证)
p3 <- records_df[records_df$spatial_score == 0 & records_df$heart_score > 0 & 
                   records_df$is_human & records_df$sc_bonus == 0 & !records_df$should_exclude, ]
p3 <- p3[order(p3$total_score, decreasing = TRUE), ]

cat("🎯 P0-人心脏空间:", nrow(p0), "条\n")
cat("🧬 P1-人心脏单细胞:", nrow(p1), "条\n")
cat("🐭 P2-小鼠心脏空间:", nrow(p2), "条\n")
cat("📊 P3-人心脏Bulk:", nrow(p3), "条\n")
cat("🚫 Exclude-非目标:", sum(records_df$should_exclude), "条\n")

# --- 保存分级结果 ---
save_grade <- function(df, filename, grade_name) {
  if (nrow(df) > 0) {
    out <- df[, c("GSE", "Title", "Summary", "Organism", "Platform", "Type", "Samples",
                  "spatial_score", "heart_score", "species_score", "total_score")]
    out$GEO_Link <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", out$GSE)
    out$人工确认 <- "待确认: 1)Supplementary文件格式 2)病理分区注释 3)样本配对信息 4)NDUFB7可检测性"
    out$推荐操作 <- ifelse(grade_name == "P0", "优先下载Series Matrix验证", 
                    ifelse(grade_name == "P1", "作为空间整合参照下载", "备选验证"))
    write.csv(out, filename, row.names = FALSE)
    cat("✅", grade_name, "已保存:", filename, "(", nrow(df), "条)\n")
    return(out)
  } else {
    cat("⚠️", grade_name, "为空，未保存文件\n")
    return(NULL)
  }
}

p0_out <- save_grade(p0, "03_results/16_spatial_from_bulk/02_P0_human_heart_spatial.csv", "P0")
p1_out <- save_grade(p1, "03_results/16_spatial_from_bulk/03_P1_human_heart_sc.csv", "P1")
p2_out <- save_grade(p2, "03_results/16_spatial_from_bulk/04_P2_mouse_heart_spatial.csv", "P2")
p3_out <- save_grade(p3, "03_results/16_spatial_from_bulk/05_P3_human_heart_bulk.csv", "P3")

# 保存全量评分
write.csv(records_df[, c("GSE", "Title", "Organism", "Type", "Samples", "spatial_score", 
                          "heart_score", "species_score", "total_score", "should_exclude")],
          "03_results/16_spatial_from_bulk/06_all_scored.csv", row.names = FALSE)

# ============================================================
# 4. 生成快速摘要报告
# ============================================================
cat("\n========== 4. 筛选摘要报告 ==========\n")

report <- file("03_results/16_spatial_from_bulk/07_screening_report.txt", "w")

writeLines("═══════════════════════════════════════════════════════════════", report)
writeLines("  GEO 2万+摘要空间转录组筛选报告", report)
writeLines(paste("  生成时间:", Sys.time()), report)
writeLines(paste("  原始记录数:", nrow(records_df)), report)
writeLines("═══════════════════════════════════════════════════════════════", report)

writeLines("\n【P0】人心脏空间转录组候选 (按总分排序)", report)
if (!is.null(p0_out) && nrow(p0_out) > 0) {
  for (i in 1:min(20, nrow(p0_out))) {
    line <- paste0(i, ". ", p0_out$GSE[i], " | 分:", p0_out$total_score[i], 
                   " | 样:", p0_out$Samples[i], " | ", substr(p0_out$Title[i], 1, 60))
    writeLines(line, report)
  }
} else {
  writeLines("  (无候选)", report)
}

writeLines("\n【P1】人心脏单细胞备选 (Top 10)", report)
if (!is.null(p1_out) && nrow(p1_out) > 0) {
  for (i in 1:min(10, nrow(p1_out))) {
    line <- paste0(i, ". ", p1_out$GSE[i], " | ", substr(p1_out$Title[i], 1, 60))
    writeLines(line, report)
  }
}

writeLines("\n【P2】小鼠心脏空间 (Top 10)", report)
if (!is.null(p2_out) && nrow(p2_out) > 0) {
  for (i in 1:min(10, nrow(p2_out))) {
    line <- paste0(i, ". ", p2_out$GSE[i], " | ", substr(p2_out$Title[i], 1, 60))
    writeLines(line, report)
  }
}

writeLines("\n【下一步行动】", report)
writeLines("1. 浏览器打开P0候选GEO链接，验证Supplementary文件", report)
writeLines("2. 确认文件格式: .h5 (SpaceRanger), .csv.gz, 或Seurat对象", report)
writeLines("3. 确认空间坐标文件存在: tissue_positions_list.csv", report)
writeLines("4. 确认病理注释: infarct/border/remote/healthy", report)
writeLines("5. 使用生成的下载脚本批量获取数据", report)
writeLines("6. 使用QC脚本验证下载完整性", report)

close(report)
cat("✅ 摘要报告已保存: 03_results/16_spatial_from_bulk/07_screening_report.txt\n")

cat("\n🎉 筛选完成！\n")
cat("【输出文件清单】\n")
cat("  01_all_parsed_records.csv      —— 全量解析结果\n")
cat("  02_P0_human_heart_spatial.csv  —— 🎯人心脏空间候选\n")
cat("  03_P1_human_heart_sc.csv       —— 🧬人心脏单细胞备选\n")
cat("  04_P2_mouse_heart_spatial.csv  —— 🐭小鼠心脏空间\n")
cat("  05_P3_human_heart_bulk.csv    —— 📊人心脏Bulk扩展\n")
cat("  06_all_scored.csv               —— 全量评分表\n")
cat("  07_screening_report.txt         —— 快速阅读报告\n")
