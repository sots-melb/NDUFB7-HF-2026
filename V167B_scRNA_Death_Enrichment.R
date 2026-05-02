#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
outdir <- file.path(PROJECT, "03_results/V167_PanDeath_Validation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V167B: GSE183852单细胞——NDUFB7-low CM中泛死亡基因富集")
message("========================================")

# 读取CM stratification
cm_file <- "03_results/V81_immediate/CM_NDUFB7_stratification.csv"
if(!file.exists(cm_file)){
  message("[FAIL] CM stratification not found: ", cm_file)
  quit(save="no", status=1)
}

cm <- fread(cm_file)
message("CM nuclei: ", nrow(cm))

# 定义各通路基因
pathway_genes <- list(
  Apoptosis = c("BCL2", "BAX", "CASP3", "CASP7", "CASP9", "BAD", "BAK1", "BID", "FAS", "TNFRSF1A"),
  Necroptosis = c("RIPK1", "RIPK3", "MLKL", "CYLD", "TNF", "TRADD", "FADD"),
  Ferroptosis_Defense = c("GPX4", "SLC7A11", "FTH1", "FTL", "NFE2L2", "SAT1"),
  Ferroptosis_Execution = c("ACSL4", "LPCAT3", "LOX", "PTGS2", "NOX1"),
  Autophagy = c("ATG5", "ATG7", "BECN1", "MAP1LC3B", "SQSTM1", "ULK1", "MTOR", "RPTOR"),
  Pyroptosis = c("NLRP3", "GSDMD", "IL1B", "IL18", "CASP1", "PYCARD")
)

# 如果metadata中有这些基因，直接比较
meta_file <- "03_results/02_tables/GSE183852_nuclei_metadata.csv"
if(file.exists(meta_file)){
  meta <- fread(meta_file)
  message("\nMetadata: ", nrow(meta), " nuclei")
  
  # 找CM
  cm_meta <- meta[meta$celltype %in% c("CM", "Cardiomyocyte", "Ventricular_CM", "Cardiomyocytes"), ]
  message("CM nuclei in metadata: ", nrow(cm_meta))
  
  if(nrow(cm_meta) > 1000 && "NDUFB7" %in% names(cm_meta)){
    # 分层
    med_nduf <- median(cm_meta$NDUFB7, na.rm = TRUE)
    cm_meta$Group <- ifelse(cm_meta$NDUFB7 < med_nduf, "NDUFB7_low", "NDUFB7_high")
    
    results <- data.frame()
    for(pw_name in names(pathway_genes)){
      genes <- pathway_genes[[pw_name]]
      available <- genes[genes %in% names(cm_meta)]
      
      if(length(available) > 0){
        message("\n[", pw_name, "] Genes found: ", paste(available, collapse = ", "))
        
        # 计算每个基因在low vs high中的差异
        for(g in available){
          low_vals <- as.numeric(cm_meta[[g]][cm_meta$Group == "NDUFB7_low"])
          high_vals <- as.numeric(cm_meta[[g]][cm_meta$Group == "NDUFB7_high"])
          low_vals <- low_vals[is.finite(low_vals)]
          high_vals <- high_vals[is.finite(high_vals)]
          
          if(length(low_vals) > 10 && length(high_vals) > 10){
            w <- wilcox.test(low_vals, high_vals, exact = FALSE)
            fc <- mean(low_vals, na.rm = TRUE) / (mean(high_vals, na.rm = TRUE) + 0.001)
            
            results <- rbind(results, data.frame(
              Pathway = pw_name,
              Gene = g,
              Low_Mean = round(mean(low_vals), 4),
              High_Mean = round(mean(high_vals), 4),
              FoldChange = round(fc, 3),
              Wilcoxon_P = format(w$p.value, digits = 2, scientific = TRUE),
              Direction = ifelse(mean(low_vals) > mean(high_vals), "UP_in_low", "DOWN_in_low"),
              stringsAsFactors = FALSE
            ))
          }
        }
      } else {
        message("\n[", pw_name, "] No genes found in metadata")
      }
    }
    
    if(nrow(results) > 0){
      fwrite(results, file.path(outdir, "V167B_scRNA_death_enrichment.csv"))
      message("\n=== 单细胞富集结果 ===")
      
      # 汇总每个通路
      summary <- results[, .(
        N_Genes = .N,
        N_UP_in_low = sum(Direction == "UP_in_low"),
        N_DOWN_in_low = sum(Direction == "DOWN_in_low"),
        Median_FC = median(FoldChange),
        Significant_Genes = sum(as.numeric(Wilcoxon_P) < 0.05, na.rm = TRUE)
      ), by = Pathway]
      print(summary[order(Median_FC, decreasing = TRUE)])
      
      # 关键判定：Ferroptosis Defense是否在low组中显著下调？
      fer_def <- results[results$Pathway == "Ferroptosis_Defense", ]
      if(nrow(fer_def) > 0){
        n_down <- sum(fer_def$Direction == "DOWN_in_low")
        n_sig_down <- sum(fer_def$Direction == "DOWN_in_low" & as.numeric(fer_def$Wilcoxon_P) < 0.05, na.rm = TRUE)
        message("\n[Ferroptosis Defense] ", n_down, "/", nrow(fer_def), " genes DOWN in NDUFB7-low CM")
        message("  Significant DOWN: ", n_sig_down, "/", nrow(fer_def))
        if(n_sig_down >= 2){
          message("  [PASS] Ferroptosis defense significantly collapsed in NDUFB7-low CM")
        }
      }
      
      # 凋亡是否也下调？
      apop <- results[results$Pathway == "Apoptosis", ]
      if(nrow(apop) > 0){
        n_up <- sum(apop$Direction == "UP_in_low")
        message("\n[Apoptosis] ", n_up, "/", nrow(apop), " genes UP in NDUFB7-low CM (execution genes should be UP)")
      }
    }
  }
} else {
  message("[FAIL] Metadata not found")
}

message("\n[DONE] V167B: ", outdir)
