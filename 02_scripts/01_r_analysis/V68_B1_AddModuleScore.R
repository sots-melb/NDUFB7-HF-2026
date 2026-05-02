suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

# --- 加载数据 ---
obj <- NULL
rds_path <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/GSE183852_Pure_CM_Subsampled.RDS"
if(file.exists(rds_path)) {
  obj <- readRDS(rds_path)
  message("▶ 从RDS加载: ", rds_path)
} else {
  robj <- "~/Projects/NDUFB7_HF_2026_04_20/01_data/01_raw_geo/GSE183852/GSE183852_DCM_Integrated.Robj"
  if(file.exists(robj)) {
    env <- new.env()
    load(robj, envir=env)
    for(v in ls(env)) {
      if(inherits(env[[v]], "Seurat")) { obj <- env[[v]]; message("▶ 从.Robj加载: ", v); break }
    }
  }
}
if(is.null(obj)) { message("❌ 无法加载GSE183852"); quit(status=1) }

message("▶ 类: ", paste(class(obj), collapse=", "))
message("▶ 维度: ", nrow(obj), " x ", ncol(obj))
message("▶ Meta.data列: ", paste(colnames(obj@meta.data), collapse=", "))

# --- 查找关键列 ---
disease_col <- NULL
for(c in c("condition","disease","group","source","cell_source","patient_group","status","orig.ident")) {
  if(c %in% colnames(obj@meta.data)) { disease_col <- c; break }
}
cell_col <- NULL
for(c in c("cell_type","celltype","CellType","cell_type_l1","cell_type_l2","seurat_clusters","predicted.id","cell_types")) {
  if(c %in% colnames(obj@meta.data)) { cell_col <- c; break }
}
if(!is.null(disease_col)) { message("▶ 疾病列: ", disease_col); print(table(obj@meta.data[[disease_col]])) }
if(!is.null(cell_col)) { message("▶ 细胞类型列: ", cell_col); print(table(obj@meta.data[[cell_col]])) }

# --- AddModuleScore: CM模块 ---
cm_genes <- c("MYH6","MYH7","TTN","ACTC1","TNNI3","RYR2","PLN","SLC8A1","MYL2","MYL3","TNNC1","TNNT2")
cm_in <- cm_genes[cm_genes %in% rownames(obj)]
message("▶ CM标记命中: ", length(cm_in), "/", length(cm_genes), " | ", paste(cm_in, collapse=", "))
if(length(cm_in)>0) {
  obj <- AddModuleScore(obj, features=list(cm_in), name="CM_Module", assay="RNA")
  message("✅ CM模块评分完成")
}

# --- AddModuleScore: 铁死亡模块 ---
ferro_genes <- c("FTL","FTH1","SLC7A11","GPX4","SAT1","ACSL4","NFE2L2","KEAP1","LPCAT3","PTGS2","ALOX15","GLS2")
ferro_in <- ferro_genes[ferro_genes %in% rownames(obj)]
message("▶ 铁死亡基因命中: ", length(ferro_in), "/", length(ferro_genes), " | ", paste(ferro_in, collapse=", "))
if(length(ferro_in)>0) {
  obj <- AddModuleScore(obj, features=list(ferro_in), name="Ferroptosis", assay="RNA")
  message("✅ 铁死亡评分完成")
}

# --- 保存 ---
out_rds <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE183852_Scored.RDS"
saveRDS(obj, out_rds)
message("✅ 保存评分对象: ", out_rds)

# --- 统计比较 ---
if(!is.null(disease_col)) {
  df <- data.frame(Disease=obj@meta.data[[disease_col]], CM_Score=obj$CM_Module1, Ferro_Score=obj$Ferroptosis1, stringsAsFactors=FALSE)
  message("▶ CM评分 by Disease:"); print(aggregate(CM_Score~Disease, data=df, FUN=function(x) c(Mean=mean(x),SD=sd(x),Median=median(x))))
  message("▶ 铁死亡评分 by Disease:"); print(aggregate(Ferro_Score~Disease, data=df, FUN=function(x) c(Mean=mean(x),SD=sd(x),Median=median(x))))
  grps <- unique(df$Disease[!is.na(df$Disease)])
  if(length(grps)==2) {
    message("▶ 双组比较:")
    print(wilcox.test(CM_Score~Disease, data=df))
    print(wilcox.test(Ferro_Score~Disease, data=df))
  }
}

# --- CM子集分析 ---
if(!is.null(cell_col)) {
  cm_labels <- grep("CM|Cardio|Myocyte|myocyte", unique(obj@meta.data[[cell_col]]), value=TRUE, ignore.case=TRUE)
  message("▶ CM标签匹配: ", paste(cm_labels, collapse=", "))
  if(length(cm_labels)>0) {
    obj_cm <- subset(obj, subset = get(cell_col) %in% cm_labels)
    message("▶ CM子集: ", ncol(obj_cm), " 细胞")
    ndufb7_cm <- tryCatch(LayerData(obj_cm, assay="RNA", layer="data")["NDUFB7", ], error=function(e) {
      tryCatch(GetAssayData(obj_cm, assay="RNA", slot="data")["NDUFB7", ], error=function(e2) NULL)
    })
    if(!is.null(ndufb7_cm)) {
      message("▶ CM中NDUFB7: Mean=", round(mean(ndufb7_cm), 4), ", %Expr=", round(mean(ndufb7_cm>0), 4))
      if("CM_Module1" %in% colnames(obj_cm@meta.data)) {
        message("▶ CM中CM评分: Mean=", round(mean(obj_cm$CM_Module1), 4))
      }
      if("Ferroptosis1" %in% colnames(obj_cm@meta.data)) {
        message("▶ CM中铁死亡评分: Mean=", round(mean(obj_cm$Ferroptosis1), 4))
      }
    }
  }
}
