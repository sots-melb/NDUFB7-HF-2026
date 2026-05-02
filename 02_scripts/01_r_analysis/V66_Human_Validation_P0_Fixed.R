suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

message("▶ 载入人类单核对象 (GSE183852)...")
rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_figures_tables/GSE183852_scored.RDS"

if(file.exists(rds_path)) {
    # 1. 内存保护：读取大文件
    full_obj <- readRDS(rds_path)
    gc()
    message(paste("✅ 读取成功。原始总细胞数:", ncol(full_obj)))
    
    # 2. 严格随机抽稀
    set.seed(42)
    keep_cells <- sample(Cells(full_obj), min(20000, ncol(full_obj)))
    obj <- subset(full_obj, cells = keep_cells)
    rm(full_obj)
    gc()
    
    # 3. 绝对提纯心肌细胞
    message("▶ 正在提纯人类心肌细胞 (TNNT2+ / MYH7+)...")
    if(all(c("TNNT2", "MYH7") %in% rownames(obj))) {
        cm_cells <- WhichCells(obj, expression = TNNT2 > 0 & MYH7 > 0)
        obj_cm <- subset(obj, cells = cm_cells)
        message(paste("✅ 提纯完成，纯净 CM 数量:", ncol(obj_cm)))
        
        # ========================================================
        # 4. 核心修复：手写极速安全打分函数，绕过 Seurat 崩溃漏洞
        # ========================================================
        safe_score <- function(seurat_obj, gene_list, score_name) {
            # 自动过滤掉矩阵中不存在的基因（如 ATP5F1A 可能被作者过滤了）
            valid_genes <- intersect(gene_list, rownames(seurat_obj))
            message(sprintf("   - 计算 %s: 命中 %d/%d 个基因", score_name, length(valid_genes), length(gene_list)))
            
            if(length(valid_genes) > 0) {
                # 提取矩阵并计算列均值 (即每个细胞的 Pathway 平均表达量)
                expr_matrix <- GetAssayData(seurat_obj, layer = "data")[valid_genes, , drop=FALSE]
                seurat_obj[[score_name]] <- colMeans(expr_matrix)
            } else {
                seurat_obj[[score_name]] <- 0
            }
            return(seurat_obj)
        }

        genes_oxphos <- c('ATP5F1A', 'COX4I1', 'CYC1', 'NDUFA4', 'SDHA')
        genes_ferro_def <- c('SLC7A11', 'GPX4', 'FTH1', 'FTL', 'NFE2L2', 'GSS')
        genes_ferro_drv <- c('ACSL4', 'PTGS2', 'ALOX15', 'TFRC')
        
        message("▶ 正在执行免疫崩溃的安全评分...")
        obj_cm <- safe_score(obj_cm, genes_oxphos, "OXPHOS_Score")
        obj_cm <- safe_score(obj_cm, genes_ferro_def, "Ferro_Def_Score")
        obj_cm <- safe_score(obj_cm, genes_ferro_drv, "Ferro_Drv_Score")
        
        # 5. 计算 NDUFB7 真实相关性
        if("NDUFB7" %in% rownames(obj_cm)) {
            expr_ndu <- as.numeric(GetAssayData(obj_cm, layer = "data")["NDUFB7", ])
            
            res_ox <- cor.test(expr_ndu, obj_cm$OXPHOS_Score)
            res_def <- cor.test(expr_ndu, obj_cm$Ferro_Def_Score)
            res_drv <- cor.test(expr_ndu, obj_cm$Ferro_Drv_Score)
            
            message("\n=======================================================")
            message("🌟 机制闭环终极结果：人类单核 (GSE183852) 纯 CM 🌟")
            message(sprintf(" 1. NDUFB7 vs OXPHOS产能:       r = %.4f (p=%.2e)", res_ox$estimate, res_ox$p.value))
            message(sprintf(" 2. NDUFB7 vs 铁死亡防线(GPX4): r = %.4f (p=%.2e)", res_def$estimate, res_def$p.value))
            message(sprintf(" 3. NDUFB7 vs 铁死亡驱动(ACSL4):r = %.4f (p=%.2e)", res_drv$estimate, res_drv$p.value))
            message("=======================================================")
        } else {
            message("❌ 严重警告：在 GSE183852 矩阵中未找到 NDUFB7 基因。")
        }
        
    } else {
        message("❌ 对象中未找到心肌标志基因。")
    }
} else {
    message("❌ 未找到原始 RDS 文件。")
}
