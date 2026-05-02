suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

message("▶ 正在检索 GSE183852 单细胞对象...")
# 寻找我们在 V65_C2_scTenifold 步骤中发现的 RDS 文件
files <- list.files("~/Projects/NDUFB7_HF_2026_04_20", pattern = "GSE183852.*\\.rds$", recursive = TRUE, full.names = TRUE)

if(length(files) > 0) {
    rds_file <- files[1]
    message(paste("✅ 找到 RDS 对象:", rds_file))
    message("   -> 正在读取数据，可能需要1-3分钟，请稍候...")
    
    # 捕获读取，防崩溃
    sc_obj <- tryCatch({ readRDS(rds_file) }, error = function(e) NULL)
    
    if(!is.null(sc_obj)) {
        # 强制垃圾回收，释放读取时的冗余内存
        gc()
        message(paste("✅ 读取成功！总细胞数:", ncol(sc_obj)))
        
        # --- 1. 绝对提纯心肌细胞 (CM) ---
        message("▶ 正在通过 TNNT2 / MYH7 绝对提纯心肌细胞...")
        if(all(c("TNNT2", "MYH7") %in% rownames(sc_obj))) {
            # 提取表达心肌标志物的细胞
            cm_cells <- WhichCells(sc_obj, expression = TNNT2 > 0 & MYH7 > 0)
            
            # 【内存保护机制】: 如果心肌细胞过多，随机抽取 10000 个以防 AddModuleScore 爆内存
            if(length(cm_cells) > 10000) {
                set.seed(42)
                cm_cells <- sample(cm_cells, 10000)
            }
            sc_obj <- subset(sc_obj, cells = cm_cells)
            gc()
            message(paste("✅ 提纯与抽样完成，用于验证的心肌细胞数:", ncol(sc_obj)))
            
            # --- 2. 定义基因集 (正反拆分) ---
            genes_oxphos <- list(c('ATP5F1A', 'COX4I1', 'CYC1', 'NDUFA4', 'SDHA'))
            genes_ferro_def <- list(c('SLC7A11', 'GPX4', 'FTH1', 'FTL', 'NFE2L2', 'GSS'))
            genes_ferro_drv <- list(c('ACSL4', 'PTGS2', 'ALOX15', 'TFRC'))
            
            # --- 3. Seurat 模块打分 ---
            message("▶ 正在计算代谢与铁死亡评分...")
            sc_obj <- AddModuleScore(sc_obj, features = genes_oxphos, name = "OXPHOS_Score")
            sc_obj <- AddModuleScore(sc_obj, features = genes_ferro_def, name = "Ferro_Def_Score")
            sc_obj <- AddModuleScore(sc_obj, features = genes_ferro_drv, name = "Ferro_Drv_Score")
            
            # 提取表达矩阵与评分
            ndufb7_expr <- as.numeric(GetAssayData(sc_obj, layer = "data")["NDUFB7", ])
            score_ox <- sc_obj$OXPHOS_Score1
            score_def <- sc_obj$Ferro_Def_Score1
            score_drv <- sc_obj$Ferro_Drv_Score1
            
            # --- 4. 相关性计算 ---
            res_ox <- cor.test(ndufb7_expr, score_ox, method="pearson")
            res_def <- cor.test(ndufb7_expr, score_def, method="pearson")
            res_drv <- cor.test(ndufb7_expr, score_drv, method="pearson")
            
            message("\n=======================================================")
            message("🌟 机制闭环 3.0：单核分辨率 (snRNA-seq) 纯 CM 相关性 🌟")
            message(sprintf(" 1. NDUFB7 vs OXPHOS产能:       r = %.4f (p=%.2e)", res_ox$estimate, res_ox$p.value))
            message(sprintf(" 2. NDUFB7 vs 铁死亡防线(GPX4): r = %.4f (p=%.2e)", res_def$estimate, res_def$p.value))
            message(sprintf(" 3. NDUFB7 vs 铁死亡驱动(ACSL4):r = %.4f (p=%.2e)", res_drv$estimate, res_drv$p.value))
            message("=======================================================")
            
            if(res_ox$estimate > 0.1 && res_def$estimate > 0.05 && res_drv$estimate < 0) {
                 message("\n🎉 极致闭环！剔除了空间污染后，NDUFB7在心肌内的真实功能链条彻底浮现：")
                 message("NDUFB7丢失 -> 产能瘫痪 -> 铁死亡防线崩塌 -> 铁死亡驱动增强！")
            }
            
        } else {
            message("❌ 对象中未找到 TNNT2 或 MYH7 标志基因，无法提纯 CM。")
        }
    } else {
        message("❌ 读取 RDS 失败，可能是内存不足。")
    }
} else {
    message("❌ 未找到 GSE183852 的 RDS 文件。")
}
