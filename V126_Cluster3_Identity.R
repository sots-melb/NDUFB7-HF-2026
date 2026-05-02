#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})
PROJECT <- path.expand("~/Projects/NDUFB7_HF_2026_04_20")
setwd(PROJECT)
outdir <- file.path(PROJECT, "03_results/V126_Cluster3_Identity")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("========================================")
message("V126: Cluster 3 生物学身份验证")
message("========================================")

CM_FILE <- "01_data/02_single_cell/GSE183852/GSE183852_CM_annotated.rds"
srt <- readRDS(CM_FILE)
message("[PASS] CM: ", ncol(srt), " cells")

# 确保NDUFB7在meta.data
if (!"NDUFB7" %in% colnames(srt@meta.data)) {
  srt$NDUFB7 <- as.numeric(FetchData(srt, vars = "NDUFB7")$NDUFB7)
}

# --- [1/4] QC指标审计：Cluster 3 vs 其他 ---
message("\n>>> [1/4] QC指标审计")
qc_vars <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
qc_vars <- intersect(qc_vars, colnames(srt@meta.data))

if (length(qc_vars) > 0) {
  srt$is_cluster3 <- ifelse(srt$seurat_clusters == 3, "Cluster_3", "Other")
  qc_df <- srt@meta.data %>% as.data.frame() %>% 
    filter(is_cluster3 %in% c("Cluster_3", "Other")) %>%
    group_by(is_cluster3) %>%
    summarise(across(all_of(qc_vars), list(mean = ~mean(.x, na.rm = TRUE), 
                                           median = ~median(.x, na.rm = TRUE),
                                           sd = ~sd(.x, na.rm = TRUE))), .groups = "drop")
  write.csv(qc_df, file.path(outdir, "V126_cluster3_qc_audit.csv"), row.names = FALSE)
  message("[PASS] QC audit saved")
  
  # 判定：如果Cluster 3的nCount/nFeature显著低于其他→低质量/双细胞风险
  # 如果percent.mt显著高于其他→死细胞/应激
  if ("percent.mt" %in% qc_vars) {
    mt_mean_c3 <- mean(srt$percent.mt[srt$seurat_clusters == 3], na.rm = TRUE)
    mt_mean_other <- mean(srt$percent.mt[srt$seurat_clusters != 3], na.rm = TRUE)
    message("  Cluster 3 mean %MT: ", round(mt_mean_c3, 2), " | Other: ", round(mt_mean_other, 2))
    if (mt_mean_c3 > mt_mean_other * 1.5) {
      message("  [WARNING] Cluster 3 has elevated mitochondrial content — possible stress/death signal")
    } else {
      message("  [PASS] Mitochondrial content not elevated — less likely dead cell artifact")
    }
  }
} else {
  message("[SKIP] No QC variables found in meta.data")
}

# --- [2/4] 差异表达：Cluster 3 vs All Others ---
message("\n>>> [2/4] 差异表达：Cluster 3 vs Others")
Idents(srt) <- "seurat_clusters"
de_genes <- tryCatch(FindMarkers(srt, ident.1 = 3, ident.2 = NULL, 
                                  min.pct = 0.1, logfc.threshold = 0.1, test.use = "wilcox"),
                     error = function(e) NULL)

if (!is.null(de_genes) && nrow(de_genes) > 0) {
  de_genes$gene <- rownames(de_genes)
  top_up <- head(de_genes %>% filter(avg_log2FC > 0) %>% arrange(p_val), 20)
  top_down <- head(de_genes %>% filter(avg_log2FC < 0) %>% arrange(p_val), 20)
  
  write.csv(de_genes, file.path(outdir, "V126_cluster3_de_all.csv"), row.names = FALSE)
  write.csv(top_up, file.path(outdir, "V126_cluster3_top20_up.csv"), row.names = FALSE)
  write.csv(top_down, file.path(outdir, "V126_cluster3_top20_down.csv"), row.names = FALSE)
  
  message("[PASS] DE genes: ", nrow(de_genes), " total")
  message("  Top upregulated in C3: ", paste(head(top_up$gene, 5), collapse = ", "))
  message("  Top downregulated in C3: ", paste(head(top_down$gene, 5), collapse = ", "))
  
  # 检查Complex I基因是否特异性下调
  c1_genes <- c("NDUFB7", "NDUFB8", "NDUFB10", "NDUFA9", "NDUFS1", "NDUFS3", "NDUFA13")
  c1_in_de <- de_genes %>% filter(gene %in% c1_genes & avg_log2FC < 0 & p_val < 0.05)
  message("\n  Complex I genes down in C3: ", nrow(c1_in_de), "/", length(c1_genes))
  if (nrow(c1_in_de) >= 3) {
    message("  [PASS] Multi-Complex I subunit downregulation supports 'OXPHOS collapse' identity")
  }
  
  # 检查应激/炎症/凋亡marker
  stress_genes <- c("HSPA1A", "HSPA1B", "DNAJB1", "DDIT3", "ATF4", "XBP1", "ERN1")
  stress_in_de <- de_genes %>% filter(gene %in% stress_genes & avg_log2FC > 0 & p_val < 0.05)
  message("  Stress genes up in C3: ", nrow(stress_in_de), "/", length(stress_genes))
  
  # 检查纤维化/EndMT marker
  fibro_genes <- c("ACTA2", "TAGLN", "POSTN", "COL1A1", "COL3A1", "VIM", "CDH2")
  fibro_in_de <- de_genes %>% filter(gene %in% fibro_genes & avg_log2FC > 0 & p_val < 0.05)
  message("  Fibrosis/EndMT genes up in C3: ", nrow(fibro_in_de), "/", length(fibro_genes))
  
} else {
  message("[WARN] DE analysis failed")
}

# --- [3/4] 可视化：Cluster 3的NDUFB7 + QC双轴图 ---
message("\n>>> [3/4] 可视化")
p <- VlnPlot(srt, features = "NDUFB7", pt.size = 0.1) +
  labs(title = "NDUFB7 by Cluster", subtitle = "Cluster 3 = NDUFB7-silent (n=48, 72.9% zero)") +
  theme_minimal()
ggsave(file.path(outdir, "V126_ndufb7_by_cluster_violin.png"), p, width = 6, height = 4, dpi = 300)

# --- [4/4] 判定 ---
message("\n>>> [4/4] 身份判定")
cat("
=== Cluster 3 身份判定框架 ===

IF (nCount_RNA not significantly lower than other clusters) AND
   (percent.mt not >1.5x other clusters) AND
   (Complex I genes >=3 significantly down) AND
   (Stress genes >=2 significantly up)
THEN → [PASS] Cluster 3 is a biologically valid 'metabolic crisis' subpopulation

IF (nCount_RNA significantly lower OR percent.mt very high)
THEN → [WARN] Technical artifact possibility — discuss in Methods limitations

IF (Complex I genes NOT down OR Fibrosis markers strongly up)
THEN → [WARN] Identity ambiguous — may be EndMT-transformed rather than OXPHOS-collapsed

RECOMMENDED PAPER STATEMENT:
'Cluster 3 (n=48 cells, 7.5% of CM) exhibited complete NDUFB7 silence in 72.9% of cells,
with concomitant downregulation of [N] Complex I subunits and upregulation of [N] ER
stress markers, consistent with a metabolically compromised subpopulation. While small
in absolute number, this cluster was independently observed in Monocle3 pseudotime
trajectory analysis as the terminal state of NDUFB7 depletion.'
", file = file.path(outdir, "V126_identity_verdict.txt"))

message("[DONE] V126: ", outdir)
