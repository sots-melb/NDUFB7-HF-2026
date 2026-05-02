suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ppcor)
})

PROJECT <- "/home/y411869/Projects/NDUFB7_HF_2026_04_20"
setwd(PROJECT)

message("========================================")
message("V171D: 四方向验证重跑")
message("========================================")

# 检查依赖
gse57338_file <- "03_results/V171_Fixed/GSE57338_gene_level_FIXED.rds"
gse59867_file <- "03_results/V171_Fixed/GSE59867_gene_level_FIXED.rds"

results <- list()

# ========================================
# D1: GSE57338 Complex I 特异性 (V168B)
# ========================================
message("")
message("=== D1: Complex I 特异性 ===")

if(file.exists(gse57338_file)) {
  gse57338 <- readRDS(gse57338_file)
  message("[PASS] GSE57338 gene matrix: ", nrow(gse57338), " x ", ncol(gse57338))
  
  # Complex I基因
  complex1_genes <- c("NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11",
                      "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS7", "NDUFS8",
                      "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA9", "NDUFA10", "NDUFA13")
  
  c1_found <- complex1_genes %in% rownames(gse57338)
  message("  Complex I found: ", sum(c1_found), "/", length(complex1_genes))
  
  if(sum(c1_found) >= 5) {
    # 泛死亡基因
    death_genes <- c("TP53", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "PARP1",
                     "MLKL", "RIPK3", "RIPK1", "PGAM5", "ACSL4", "GPX4", "SLC7A11",
                     "GSDMD", "GSDME", "IL1B", "NLRP3", "AIM2")
    death_found <- death_genes %in% rownames(gse57338)
    message("  Death genes found: ", sum(death_found), "/", length(death_genes))
    
    if(sum(death_found) >= 5) {
      # Spearman相关
      c1_mat <- gse57338[complex1_genes[c1_found], , drop=FALSE]
      death_mat <- gse57338[death_genes[death_found], , drop=FALSE]
      
      # NDUFB7与每个死亡基因的相关
      ndufb7_vec <- gse57338["NDUFB7", ]
      cor_results <- data.frame(
        DeathGene = death_genes[death_found],
        Spearman_r = NA,
        p_value = NA
      )
      
      for(i in 1:nrow(death_mat)) {
        g <- death_genes[death_found][i]
        test <- cor.test(ndufb7_vec, death_mat[g, ], method="spearman", exact=FALSE)
        cor_results$Spearman_r[i] <- test$estimate
        cor_results$p_value[i] <- test$p.value
      }
      
      # 排序：NDUFB7是否为Complex I中与死亡最相关的？
      c1_cor <- data.frame(
        ComplexIGene = complex1_genes[c1_found],
        MeanAbsCor = NA
      )
      for(i in 1:nrow(c1_mat)) {
        gene <- complex1_genes[c1_found][i]
        cors <- sapply(1:nrow(death_mat), function(j) {
          abs(cor.test(c1_mat[gene,], death_mat[j,], method="spearman", exact=FALSE)$estimate)
        })
        c1_cor$MeanAbsCor[i] <- mean(cors)
      }
      
      c1_cor <- c1_cor[order(-c1_cor$MeanAbsCor),]
      ndufb7_rank <- which(c1_cor$ComplexIGene == "NDUFB7")
      
      message("  NDUFB7 rank in Complex I (death correlation): ", ndufb7_rank, "/", nrow(c1_cor))
      message("  Top 3: ", paste(head(c1_cor$ComplexIGene, 3), collapse=", "))
      
      write.csv(cor_results, "03_results/V171_Fixed/D1_Complex1_Death_Cor.csv", row.names=FALSE)
      write.csv(c1_cor, "03_results/V171_Fixed/D1_Complex1_Ranking.csv", row.names=FALSE)
      
      results$D1 <- list(status="PASS", ndufb7_rank=ndufb7_rank, total=nrow(c1_cor))
    } else {
      message("[SKIP] Too few death genes")
      results$D1 <- list(status="SKIP", reason="Too few death genes")
    }
  } else {
    message("[SKIP] Too few Complex I genes")
    results$D1 <- list(status="SKIP", reason="Too few Complex I genes")
  }
} else {
  message("[SKIP] GSE57338 gene matrix not available")
  results$D1 <- list(status="SKIP", reason="No GSE57338 matrix")
}

# ========================================
# D2: 偏相关控制ROS (V168D)
# ========================================
message("")
message("=== D2: 偏相关 (控制ROS) ===")

if(file.exists(gse57338_file)) {
  gse57338 <- readRDS(gse57338_file)
  
  ros_genes <- c("NOX1", "NOX2", "NOX4", "SOD1", "SOD2", "CAT", "PRDX1", "TXNRD1", "GCLC", "GCLM")
  death_genes2 <- c("ACSL4", "GPX4", "SLC7A11", "TP53", "BAX", "CASP3", "MLKL", "RIPK3", "GSDMD")
  
  ros_found <- ros_genes %in% rownames(gse57338)
  death_found2 <- death_genes2 %in% rownames(gse57338)
  
  message("  ROS genes found: ", sum(ros_found), "/", length(ros_genes))
  message("  Death genes found: ", sum(death_found2), "/", length(death_genes2))
  
  if("NDUFB7" %in% rownames(gse57338) && sum(ros_found) >= 3 && sum(death_found2) >= 3) {
    # 构建数据框
    all_genes <- c("NDUFB7", ros_genes[ros_found], death_genes2[death_found2])
    df <- as.data.frame(t(gse57338[all_genes, ]))
    
    # 简单偏相关：NDUFB7 vs Death controlling ROS
    # 使用ppcor::pcor.test
    pcor_results <- data.frame(
      DeathGene = character(),
      Simple_r = numeric(),
      Partial_r = numeric(),
      p_value = numeric(),
      Attenuation = numeric()  # (simple - partial)/simple
    )
    
    for(d in death_genes2[death_found2]) {
      if(d == "NDUFB7") next
      
      # Simple correlation
      simple <- cor.test(df$NDUFB7, df[[d]], method="spearman", exact=FALSE)
      
      # Partial correlation controlling all ROS genes
      ros_vars <- ros_genes[ros_found]
      if(length(ros_vars) >= 1) {
        tryCatch({
          pc <- pcor.test(df$NDUFB7, df[[d]], df[, ros_vars, drop=FALSE], method="spearman")
          att <- (abs(simple$estimate) - abs(pc$estimate)) / abs(simple$estimate)
          
          pcor_results <- rbind(pcor_results, data.frame(
            DeathGene = d,
            Simple_r = simple$estimate,
            Partial_r = pc$estimate,
            p_value = pc$p.value,
            Attenuation = att
          ))
        }, error = function(e) {
          message("  [WARN] pcor failed for ", d, ": ", conditionMessage(e))
        })
      }
    }
    
    if(nrow(pcor_results) > 0) {
      message("  Partial correlation results:")
      print(pcor_results)
      
      # 判定：Attenuation > 0.3 表示ROS介导
      mediated <- sum(pcor_results$Attenuation > 0.3, na.rm=TRUE)
      message("  ROS-mediated (Attenuation>0.3): ", mediated, "/", nrow(pcor_results))
      
      write.csv(pcor_results, "03_results/V171_Fixed/D2_Partial_Correlation.csv", row.names=FALSE)
      results$D2 <- list(status="PASS", mediated=mediated, total=nrow(pcor_results))
    } else {
      message("[SKIP] No partial correlation results")
      results$D2 <- list(status="SKIP", reason="No results")
    }
  } else {
    message("[SKIP] Missing genes for partial correlation")
    results$D2 <- list(status="SKIP", reason="Missing genes")
  }
} else {
  message("[SKIP] No GSE57338 matrix")
  results$D2 <- list(status="SKIP", reason="No matrix")
}

# ========================================
# D3: GSE59867 ACSL4/GPX4比值 (V168C)
# ========================================
message("")
message("=== D3: ACSL4/GPX4比值验证 ===")

if(file.exists(gse59867_file) && file.exists(gse57338_file)) {
  gse59867 <- readRDS(gse59867_file)
  gse57338 <- readRDS(gse57338_file)
  
  if(all(c("ACSL4", "GPX4") %in% rownames(gse59867))) {
    # 计算比值
    ratio59867 <- gse59867["ACSL4", ] / (gse59867["GPX4", ] + 0.001)
    
    # 与GSE57338的NDUFB7做相关（跨队列）
    # 注意：这是队列级相关，需要谨慎解释
    message("  GSE59867 ACSL4/GPX4 ratio calculated")
    message("  Range: ", round(min(ratio59867), 3), " - ", round(max(ratio59867), 3))
    
    # 保存
    write.csv(data.frame(Sample=names(ratio59867), ACSL4_GPX4_Ratio=ratio59867),
              "03_results/V171_Fixed/D3_ACSL4_GPX4_Ratio.csv", row.names=FALSE)
    results$D3 <- list(status="PASS")
  } else {
    message("[SKIP] ACSL4 or GPX4 not found in GSE59867")
    results$D3 <- list(status="SKIP", reason="Missing ACSL4/GPX4")
  }
} else {
  message("[SKIP] Matrices not available")
  results$D3 <- list(status="SKIP", reason="No matrices")
}

# ========================================
# D4: GSE157282 DOX模型 (占位)
# ========================================
message("")
message("=== D4: GSE157282 DOX模型 ===")
message("  [PENDING] Waiting for GSE157282 matrix construction")

results$D4 <- list(status="PENDING")

# ========================================
# 汇总
# ========================================
message("")
message("========================================")
message("V171D 四方向验证汇总")
message("========================================")
for(name in names(results)) {
  r <- results[[name]]
  message(sprintf("  %s: %s", name, r$status))
}

message("[DONE] All results in 03_results/V171_Fixed/")
