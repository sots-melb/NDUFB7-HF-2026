suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

# 读取已冻结的纯心肌细胞对象
rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS"
if(file.exists(rds_path)) {
  obj <- readRDS(rds_path)
  
  # 如果已有 Condition 列，寻找 DCM 和 Donor (或者类似的分组)
  if("Condition" %in% colnames(obj@meta.data)) {
    Idents(obj) <- "Condition"
    message("▶ 正在计算 DCM vs Donor 的差异基因 (这可能需要几分钟)...")
    
    # 获取表达上调和下调的靶点
    degs <- FindMarkers(obj, ident.1 = "DCM", ident.2 = "Donor", logfc.threshold = 0.25, min.pct = 0.1)
    degs$Gene <- rownames(degs)
    
    # 过滤出显著基因 (p_val_adj < 0.05)
    sig_degs <- degs %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
    
    if(nrow(sig_degs) > 0) {
      up_genes <- sig_degs %>% filter(avg_log2FC > 0) %>% head(150) %>% pull(Gene)
      down_genes <- sig_degs %>% filter(avg_log2FC < 0) %>% head(150) %>% pull(Gene)
      
      # 分别保存为 cMap 接受的格式 (.txt)
      up_file <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/cMap_Input/DCM_CM_Up_Genes.txt"
      down_file <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/cMap_Input/DCM_CM_Down_Genes.txt"
      
      write.table(up_genes, up_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
      write.table(down_genes, down_file, row.names=FALSE, col.names=FALSE, quote=FALSE)
      
      message("✅ 成功提取 cMap 签名！")
      message("   上调基因: ", length(up_genes), " 个 (已保存至 ", up_file, ")")
      message("   下调基因: ", length(down_genes), " 个 (已保存至 ", down_file, ")")
      message("💡 下一步操作: 请将这两个 txt 文件的内容复制粘贴到 https://clue.io/query 的输入框中进行药物筛选！")
    } else {
      message("⚠️ 未找到显著差异基因。")
    }
  } else {
    message("⚠️ Seurat 对象中未找到 Condition 分组列。")
  }
} else {
  message("❌ 找不到 GSE183852_Pure_CM_Subsampled.RDS。请确认路径。")
}
