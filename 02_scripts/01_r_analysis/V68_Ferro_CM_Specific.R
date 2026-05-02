suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

obj <- readRDS("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE183852_Scored.RDS")
message("▶ 总细胞: ", ncol(obj))

# 策略1: 基于seurat_clusters识别CM
clusters <- table(obj$seurat_clusters)
message("▶ 全部Cluster大小:")
print(clusters)

# 查找CM相关cluster（基于CM_score1均值）
cm_scores <- aggregate(CM_score1~seurat_clusters, data=obj@meta.data, FUN=mean)
cm_scores <- cm_scores[order(-cm_scores$CM_score1), ]
message("▶ Cluster按CM_score排序:")
print(head(cm_scores, 10))

# 选择CM_score最高的cluster作为CM cluster
top_cm_clusters <- head(cm_scores$seurat_clusters, 3)
message("▶ Top 3 CM clusters: ", paste(top_cm_clusters, collapse=", "))

# 提取CM细胞
cm_cells <- colnames(obj)[obj$seurat_clusters %in% top_cm_clusters]
message("▶ CM细胞数: ", length(cm_cells))

if(length(cm_cells)>50) {
  obj_cm <- subset(obj, cells=cm_cells)
  df_cm <- data.frame(
    Condition=obj_cm$Condition,
    NDUFB7=obj_cm$NDUFB7_expr,
    Ferro=obj_cm$Ferro_score1,
    OXPHOS=obj_cm$OXPHOS_Score,
    CM_Score=obj_cm$CM_score1,
    stringsAsFactors=FALSE
  )
  
  message("▶ CM子集统计:")
  print(aggregate(.~Condition, data=df_cm, FUN=function(x) c(Mean=mean(x),SD=sd(x),Median=median(x))))
  
  if(length(unique(df_cm$Condition))==2) {
    message("▶ CM子集Wilcoxon:")
    print(wilcox.test(NDUFB7~Condition, data=df_cm))
    print(wilcox.test(Ferro~Condition, data=df_cm))
    print(wilcox.test(OXPHOS~Condition, data=df_cm))
    print(wilcox.test(CM_Score~Condition, data=df_cm))
  }
  
  # CM内部NDUFB7与铁死亡的相关性
  message("▶ CM内部NDUFB7 vs Ferro相关: r=", round(cor(df_cm$NDUFB7, df_cm$Ferro, use="complete.obs"), 3))
  message("▶ CM内部NDUFB7 vs OXPHOS相关: r=", round(cor(df_cm$NDUFB7, df_cm$OXPHOS, use="complete.obs"), 3))
  
  # NDUFB7分层分析
  df_cm$NDUFB7_Level <- ifelse(df_cm$NDUFB7==0, "Zero",
                        ifelse(df_cm$NDUFB7<=median(df_cm$NDUFB7[df_cm$NDUFB7>0]), "Low", "High"))
  message("▶ CM中NDUFB7分层 vs Ferro:")
  print(aggregate(Ferro~NDUFB7_Level+Condition, data=df_cm, FUN=function(x) c(Mean=mean(x),SD=sd(x))))
}

# 策略2: 非CM细胞中的铁死亡
non_cm_cells <- colnames(obj)[!obj$seurat_clusters %in% top_cm_clusters]
message("▶ 非CM细胞数: ", length(non_cm_cells))
if(length(non_cm_cells)>50) {
  obj_non <- subset(obj, cells=non_cm_cells)
  df_non <- data.frame(Condition=obj_non$Condition, Ferro=obj_non$Ferro_score1, NDUFB7=obj_non$NDUFB7_expr, stringsAsFactors=FALSE)
  message("▶ 非CM子集Ferro统计:")
  print(aggregate(Ferro~Condition, data=df_non, FUN=function(x) c(Mean=mean(x),SD=sd(x),Median=median(x))))
  if(length(unique(df_non$Condition))==2) {
    message("▶ 非CM Ferro Wilcoxon:")
    print(wilcox.test(Ferro~Condition, data=df_non))
  }
}

# 策略3: 纤维化相关cluster（Ferro_Def最高）
ferro_clusters <- aggregate(Ferro_Def_Score~seurat_clusters, data=obj@meta.data, FUN=mean)
ferro_clusters <- ferro_clusters[order(-ferro_clusters$Ferro_Def_Score), ]
message("▶ Cluster按Ferro_Def排序:")
print(head(ferro_clusters, 10))
