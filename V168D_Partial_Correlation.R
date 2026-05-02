#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ppcor) })
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V168_PanDeath_Solidity/V168D_Partial_Cor")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V168D: NDUFB7-死亡通路偏相关（控制ROS）")
message("========================================")

# === 1. 读取表达矩阵 ===
expr_file <- "03_results/02_tables/GSE57338_expression_matrix_annotated.csv"
if(!file.exists(expr_file)){
  expr_file <- "03_results/02_tables/GSE57338_expression_matrix_raw.csv"
}

if(!file.exists(expr_file)){
  message("[FAIL] Expression matrix not found")
  quit(save="no", status=1)
}

message("\n[1/3] 读取表达矩阵...")
df <- fread(expr_file, header = TRUE)
gene_col <- names(df)[1]
is_gene_row <- any(grepl("NDUFB7|ENSG|Gene", df[[1]][1:5], ignore.case = TRUE))

if(!is_gene_row){
  message("[FAIL] 需要行=基因的矩阵格式")
  quit(save="no", status=1)
}

# === 2. 提取目标基因 ===
target_genes <- c("NDUFB7", "NOX4", "SOD2", "CAT", "GPX1", "DUOX2")
death_genes <- list(
  Apoptosis = c("BCL2", "BAX", "CASP3", "CASP7", "CASP9"),
  Necroptosis = c("RIPK1", "RIPK3", "MLKL"),
  Ferroptosis_Defense = c("GPX4", "SLC7A11", "FTH1", "FTL"),
  Ferroptosis_Execution = c("ACSL4", "LPCAT3", "LOX")
)

all_needed <- unique(c(target_genes, unlist(death_genes)))
found <- all_needed[tolower(all_needed) %in% tolower(df[[gene_col]])]
message("  目标基因找到: ", length(found), "/", length(all_needed))
message("  缺失: ", paste(setdiff(all_needed, found), collapse = ", "))

if(length(found) < 10){
  message("[FAIL] Too few genes found for partial correlation")
  quit(save="no", status=1)
}

# 提取表达值
num_cols <- setdiff(names(df), gene_col)
idx <- match(tolower(found), tolower(df[[gene_col]]))
expr_mat <- as.matrix(df[idx, ..num_cols])
rownames(expr_mat) <- found
expr_mat <- t(apply(expr_mat, 2, as.numeric))  # 基因×样本
colnames(expr_mat) <- num_cols

# === 3. 计算样本级通路评分 ===
message("\n[2/3] 计算样本级通路评分...")
pw_scores <- data.frame(Sample = num_cols, stringsAsFactors = FALSE)
for(pw in names(death_genes)){
  g <- death_genes[[pw]][tolower(death_genes[[pw]]) %in% tolower(found)]
  if(length(g) > 0){
    vals <- expr_mat[tolower(g), , drop = FALSE]
    if(is.matrix(vals) && nrow(vals) > 0){
      pw_scores[[pw]] <- colMeans(vals, na.rm = TRUE)
      message("  ", pw, ": ", length(g), " genes, mean score range: ", 
              round(min(pw_scores[[pw]], na.rm = TRUE), 2), " - ", 
              round(max(pw_scores[[pw]], na.rm = TRUE), 2))
    }
  }
}

# === 4. 偏相关分析 ===
message("\n[3/3] 偏相关分析（控制ROS基因）...")
ros_genes <- c("NOX4", "SOD2", "CAT", "GPX1")[c("NOX4", "SOD2", "CAT", "GPX1") %in% found]

if(length(ros_genes) >= 2){
  # 构建数据框
  cor_df <- data.frame(
    NDUFB7 = as.numeric(expr_mat["NDUFB7", ]),
    stringsAsFactors = FALSE
  )
  for(pw in setdiff(names(pw_scores), "Sample")){
    cor_df[[pw]] <- as.numeric(pw_scores[[pw]])
  }
  for(rg in ros_genes){
    cor_df[[rg]] <- as.numeric(expr_mat[rg, ])
  }
  
  # 移除NA
  cor_df <- cor_df[complete.cases(cor_df), ]
  message("  Complete cases: ", nrow(cor_df))
  
  if(nrow(cor_df) > 50){
    results <- data.frame()
    for(pw in setdiff(names(pw_scores), "Sample")){
      if(pw %in% names(cor_df)){
        # 简单偏相关：控制所有ROS基因
        tryCatch({
          pc <- pcor.test(cor_df$NDUFB7, cor_df[[pw]], cor_df[, ros_genes, drop = FALSE])
          results <- rbind(results, data.frame(
            Pathway = pw,
            ROS_Control = paste(ros_genes, collapse = "+"),
            Partial_R = round(pc$estimate, 3),
            Partial_P = format(pc$p.value, digits = 2, scientific = TRUE),
            N = nrow(cor_df),
            stringsAsFactors = FALSE
          ))
        }, error = function(e) message("  ", pw, " failed: ", conditionMessage(e)))
      }
    }
    
    if(nrow(results) > 0){
      fwrite(results, file.path(outdir, "V168D_partial_correlation.csv"))
      message("\n=== 偏相关结果（控制ROS）===")
      print(results[order(abs(Partial_R), decreasing = TRUE), ])
      
      # 与简单相关比较
      simple <- data.frame()
      for(pw in setdiff(names(pw_scores), "Sample")){
        if(pw %in% names(cor_df)){
          r <- cor(cor_df$NDUFB7, cor_df[[pw]], method = "spearman")
          simple <- rbind(simple, data.frame(Pathway = pw, Simple_Rho = round(r, 3), stringsAsFactors = FALSE))
        }
      }
      
      comp <- merge(results[, c("Pathway", "Partial_R")], simple, by = "Pathway")
      comp$Attenuation <- round((1 - abs(comp$Partial_R) / abs(comp$Simple_Rho)) * 100, 1)
      message("\n=== 简单相关 vs 偏相关 ===")
      print(comp[order(comp$Attenuation, decreasing = TRUE), ])
      
      # 判定：如果控制ROS后相关大幅衰减→说明NDUFB7通过ROS介导死亡
      high_atten <- comp$Attenuation[comp$Attenuation > 50 & !is.na(comp$Attenuation)]
      if(length(high_atten) > 0){
        message("\n[WARN] ", length(high_atten), " pathways show >50% attenuation after ROS control")
        message("  → NDUFB7 may mediate pan-death primarily through ROS")
      } else {
        message("\n[PASS] NDUFB7-death associations are largely independent of ROS genes")
      }
    }
  }
} else {
  message("[WARN] Too few ROS genes found (", length(ros_genes), ") — need NOX4/SOD2/CAT/GPX1 for partial correlation")
  message("  标记为Revision期任务")
}

message("\n[DONE] V168D: ", outdir)
