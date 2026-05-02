suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

message("▶ 载入人类单核对象 (GSE183852)...")
rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_figures_tables/GSE183852_scored.RDS"

if(file.exists(rds_path)) {
    # 1. 内存保护：读取大文件
    full_obj <- readRDS(rds_path)
    gc()
    message(paste("✅ 读取成功。原始总细胞数:", ncol(full_obj)))
    
    # 2. 严格随机抽稀至 20,000 细胞以避开 OOM
    set.seed(42)
    keep_cells <- sample(Cells(full_obj), min(20000, ncol(full_obj)))
    obj <- subset(full_obj, cells = keep_cells)
    rm(full_obj) # 立即释放大内存
    gc()
    message(paste("✅ 降采样完成，当前分析样本量:", ncol(obj)))
    
    # 3. 绝对提纯心肌细胞 (CM)
    message("▶ 正在提纯人类心肌细胞 (TNNT2+ / MYH7+)...")
    if(all(c("TNNT2", "MYH7") %in% rownames(obj))) {
        cm_cells <- WhichCells(obj, expression = TNNT2 > 0 & MYH7 > 0)
        obj_cm <- subset(obj, cells = cm_cells)
        message(paste("✅ 提纯完成，纯净 CM 数量:", ncol(obj_cm)))
        
        # 4. 定义 V65 确立的机制基因集 (正反拆分)
        genes_oxphos <- list(c('ATP5F1A', 'COX4I1', 'CYC1', 'NDUFA4', 'SDHA'))
        genes_ferro_def <- list(c('SLC7A11', 'GPX4', 'FTH1', 'FTL', 'NFE2L2', 'GSS'))
        genes_ferro_drv <- list(c('ACSL4', 'PTGS2', 'ALOX15', 'TFRC'))
        
        # 5. 运行模块打分
        message("▶ 正在计算代谢与铁死亡评分...")
        obj_cm <- AddModuleScore(obj_cm, features = genes_oxphos, name = "OXPHOS_Score")
        obj_cm <- AddModuleScore(obj_cm, features = genes_ferro_def, name = "Ferro_Def_Score")
        obj_cm <- AddModuleScore(obj_cm, features = genes_ferro_drv, name = "Ferro_Drv_Score")
        
        # 6. 计算 NDUFB7 真实相关性
        expr_ndu <- as.numeric(GetAssayData(obj_cm, layer = "data")["NDUFB7", ])
        res_ox <- cor.test(expr_ndu, obj_cm$OXPHOS_Score1)
        res_def <- cor.test(expr_ndu, obj_cm$Ferro_Def_Score1)
        res_drv <- cor.test(expr_ndu, obj_cm$Ferro_Drv_Score1)
        
        # 7. 打印终极大满贯结果
        message("\n=======================================================")
        message("🌟 机制闭环最终结果：人类单核 (snRNA-seq) 纯 CM 🌟")
        message(sprintf(" 1. NDUFB7 vs OXPHOS产能:       r = %.4f (p=%.2e)", res_ox$estimate, res_ox$p.value))
        message(sprintf(" 2. NDUFB7 vs 铁死亡防线(GPX4): r = %.4f (p=%.2e)", res_def$estimate, res_def$p.value))
        message(sprintf(" 3. NDUFB7 vs 铁死亡驱动(ACSL4):r = %.4f (p=%.2e)", res_drv$estimate, res_drv$p.value))
        message("=======================================================")
        
        # 保存轻量化对象供 Figure 使用
        saveRDS(obj_cm, "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS")
        message("🎉 轻量化验证对象已保存。")
        
    } else {
        message("❌ 对象中未找到标志基因。")
    }
} else {
    message("❌ 未找到原始 RDS 文件。")
}
