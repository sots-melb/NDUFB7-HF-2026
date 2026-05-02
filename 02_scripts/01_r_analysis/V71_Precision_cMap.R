suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS"
if(file.exists(rds_path)) {
  obj <- readRDS(rds_path)
  
  # 提取 NDUFB7 表达量
  ndufb7_expr <- tryCatch({
    GetAssayData(obj, slot="data")["NDUFB7", ]
  }, error = function(e) {
    LayerData(obj, assay="RNA", layer="data")["NDUFB7", ]
  })
  
  # 根据 NDUFB7 的表达量将纯心肌细胞分为高/低表达组
  q25 <- quantile(ndufb7_expr[ndufb7_expr > 0], 0.25)
  q75 <- quantile(ndufb7_expr[ndufb7_expr > 0], 0.75)
  
  obj$NDUFB7_Group <- ifelse(ndufb7_expr <= q25, "Low", 
                      ifelse(ndufb7_expr >= q75, "High", "Mid"))
                      
  # 过滤掉中间的细胞，只比较极端情况
  obj_extreme <- subset(obj, subset = NDUFB7_Group %in% c("Low", "High"))
  Idents(obj_extreme) <- "NDUFB7_Group"
  
  message("▶ 成功划分 NDUFB7 分组：Low (n=", sum(obj$NDUFB7_Group=="Low"), ") vs High (n=", sum(obj$NDUFB7_Group=="High"), ")")
  message("▶ 正在计算特异性差异基因 (Low vs High)...")
  
  # 放宽阈值，抓取单细胞中的微弱但一致的信号
  degs <- FindMarkers(obj_extreme, ident.1 = "Low", ident.2 = "High", logfc.threshold = 0.1, min.pct = 0.05)
  degs$Gene <- rownames(degs)
  
  # 过滤出显著基因
  sig_degs <- degs %>% filter(p_val_adj < 0.05)
  
  if(nrow(sig_degs) > 0) {
    # 按照 FoldChange 排序，提取最靠前的 150 个
    up_genes <- sig_degs %>% filter(avg_log2FC > 0) %>% arrange(desc(avg_log2FC)) %>% head(150) %>% pull(Gene)
    down_genes <- sig_degs %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC) %>% head(150) %>% pull(Gene)
    
    up_file <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/cMap_Input/NDUFB7_Specific_Up_Genes.txt"
    down_file <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/cMap_Input/NDUFB7_Specific_Down_Genes.txt"
    
    write.table(up_genes, up_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(down_genes, down_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    message("✅ 完美提取 NDUFB7 专属靶向签名！")
    message("   上调基因: ", length(up_genes), " 个 (已保存至 ", up_file, ")")
    message("   下调基因: ", length(down_genes), " 个 (已保存至 ", down_file, ")")
  } else {
    message("⚠️ 未找到显著差异基因。")
  }
}
