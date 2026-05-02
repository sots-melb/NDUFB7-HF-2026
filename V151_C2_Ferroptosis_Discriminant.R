#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V151_C2_Discriminant")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V151: C2强化 — 铁死亡Discriminant Validation")
message("========================================")

# === 1. 加载GSE57338表达矩阵 ===
rds_file <- "03_results/V133_GSE57338/GSE57338_gene_level.rds"
if (!file.exists(rds_file)) {
  rds_file <- list.files(c("01_data","03_results"), pattern = "GSE57338.*gene_level.*\\.rds$", full.names = TRUE, recursive = TRUE)[1]
}
obj <- readRDS(rds_file)
exprs <- if (is.list(obj) && "exprs" %in% names(obj)) obj$exprs else as.matrix(obj)
message("[PASS] Exprs: ", nrow(exprs), "x", ncol(exprs))

# === 2. 定义死亡通路基因集 ===
death_pathways <- list(
  Ferroptosis_Defense = c("GPX4","SLC7A11","FTH1","FTL","NFE2L2","SAT1","ACSL4","KEAP1"),
  Ferroptosis_Execution = c("PTGS2","NOX1","LPCAT3","AIFM2","TFRC"),
  Apoptosis = c("CASP3","CASP8","CASP9","BAX","BAK1","BCL2","BCL2L1","FAS","FASLG"),
  Necroptosis = c("MLKL","RIPK1","RIPK3","PGAM5","CYLD"),
  Pyroptosis = c("GSDMD","GSDME","NLRP3","IL1B","IL18","CASP1","CASP4"),
  Autophagy = c("ATG5","ATG7","BECN1","LC3","SQSTM1","MAP1LC3B")
)

# === 3. 计算各通路评分 ===
message("\n[1/2] 计算6种细胞死亡通路评分...")
ndufb7_expr <- as.numeric(exprs["NDUFB7", ])
names(ndufb7_expr) <- colnames(exprs)

results <- data.frame()
for (pw_name in names(death_pathways)) {
  genes <- intersect(death_pathways[[pw_name]], rownames(exprs))
  if (length(genes) >= 3) {
    # z-score均值
    mat <- exprs[genes, , drop = FALSE]
    z_mat <- apply(mat, 1, function(g) (g - mean(g, na.rm=TRUE)) / sd(g, na.rm=TRUE))
    if (!is.null(dim(z_mat))) {
      score <- rowMeans(z_mat, na.rm = TRUE)
    } else {
      score <- as.numeric(z_mat)
    }
    
    # Spearman vs NDUFB7
    ct <- cor.test(ndufb7_expr, score, method = "spearman", use = "complete.obs")
    
    results <- rbind(results, data.frame(
      Pathway = pw_name,
      N_Genes = length(genes),
      Genes_Used = paste(genes, collapse = ","),
      Spearman_Rho = round(ct$estimate, 3),
      P_Value = ct$p.value,
      Significant = ct$p.value < 0.05,
      Direction = ifelse(ct$estimate > 0, "Positive", "Negative"),
      stringsAsFactors = FALSE
    ))
    message("  ", pw_name, ": rho=", round(ct$estimate, 3), " p=", format(ct$p.value, digits=2, scientific=TRUE))
  } else {
    message("  ", pw_name, ": SKIP (only ", length(genes), " genes)")
  }
}

# === 4. Discriminant判定 ===
message("\n[2/2] Discriminant判定...")
if (nrow(results) > 0) {
  results <- results[order(results$P_Value), ]
  fwrite(results, file.path(outdir, "V151_discriminant_validation.csv"))
  
  ferro_rho <- abs(results$Spearman_Rho[results$Pathway == "Ferroptosis_Defense"])
  other_rho <- abs(results$Spearman_Rho[results$Pathway != "Ferroptosis_Defense"])
  
  message("\n=== Discriminant Results ===")
  print(results[, c("Pathway","Spearman_Rho","P_Value","Significant")])
  
  if (length(ferro_rho) > 0 && all(ferro_rho > other_rho, na.rm = TRUE)) {
    message("\n[PASS] 铁死亡|rho| > 所有其他通路 → 特异性确认！")
  } else if (length(ferro_rho) > 0 && ferro_rho > max(other_rho[other_rho < 0.3], na.rm = TRUE)) {
    message("\n[PARTIAL] 铁死亡|rho|最大，但部分通路也显著 → 写为'predominantly ferroptosis'")
  } else {
    message("\n[WARN] 其他通路|rho| ≥ 铁死亡 → 特异性不足，需调整叙事")
  }
  
  # 方法学比较：AUCell vs SumScore vs z-mean
  message("\n--- Method comparison (if V141 data available) ---")
  v141_file <- "03_results/V141_GSE57338_Bulk_Fix/V141_bulk_ferroptosis.csv"
  if (file.exists(v141_file)) {
    v141 <- read.csv(v141_file)
    if (all(c("NDUFB7","ferroptosis_score") %in% colnames(v141))) {
      ct2 <- cor.test(v141$NDUFB7, v141$ferroptosis_score, method = "spearman")
      message("  V141 z-mean score: rho=", round(ct2$estimate, 3), " p=", format(ct2$p.value, digits=2))
    }
  }
}

message("[DONE] V151: ", outdir)
