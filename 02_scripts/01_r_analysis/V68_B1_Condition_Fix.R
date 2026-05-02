suppressMessages(library(Seurat))
obj <- readRDS("~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE183852_Scored.RDS")
message("▶ 维度: ", nrow(obj), " x ", ncol(obj))
message("▶ Condition分布:")
print(table(obj$Condition, useNA="ifany"))
message("▶ seurat_clusters分布:")
print(table(obj$seurat_clusters, useNA="ifany"))

# 按Condition比较
df <- data.frame(Condition=obj$Condition, CM_Score=obj$CM_score1, Ferro_Score=obj$Ferro_score1, OXPHOS=obj$OXPHOS_Score, Ferro_Def=obj$Ferro_Def_Score, NDUFB7=obj$NDUFB7_expr, stringsAsFactors=FALSE)
message("▶ 各组统计:")
for(col in c("CM_Score","Ferro_Score","OXPHOS","Ferro_Def","NDUFB7")) {
  message("  --- ", col, " ---")
  print(aggregate(as.formula(paste(col,"~Condition")), data=df, FUN=function(x) c(Mean=mean(x),SD=sd(x),Median=median(x),N=length(x))))
}
grps <- unique(df$Condition[!is.na(df$Condition)])
if(length(grps)==2) {
  message("▶ 双组Wilcoxon:")
  for(col in c("CM_Score","Ferro_Score","OXPHOS","Ferro_Def","NDUFB7")) {
    f <- as.formula(paste(col,"~Condition"))
    res <- wilcox.test(f, data=df)
    message(sprintf("  %s: W=%.1f, p=%.4f", col, res$statistic, res$p.value))
  }
}

# CM子集（基于CM_score1 top分位数）
cm_cutoff <- quantile(df$CM_Score, 0.8, na.rm=TRUE)
cm_cells <- colnames(obj)[obj$CM_score1 >= cm_cutoff]
message("▶ CM-like细胞数 (CM_score top 20%): ", length(cm_cells))
if(length(cm_cells)>10) {
  obj_cm <- subset(obj, cells=cm_cells)
  df_cm <- data.frame(Condition=obj_cm$Condition, NDUFB7=obj_cm$NDUFB7_expr, Ferro=obj_cm$Ferro_score1, OXPHOS=obj_cm$OXPHOS_Score, stringsAsFactors=FALSE)
  message("▶ CM-like子集统计:")
  print(aggregate(NDUFB7~Condition, data=df_cm, FUN=function(x) c(Mean=mean(x),SD=sd(x),Median=median(x))))
  print(aggregate(Ferro~Condition, data=df_cm, FUN=function(x) c(Mean=mean(x),SD=sd(x),Median=median(x))))
  print(aggregate(OXPHOS~Condition, data=df_cm, FUN=function(x) c(Mean=mean(x),SD=sd(x),Median=median(x))))
  if(length(unique(df_cm$Condition))==2) {
    message("▶ CM-like Wilcoxon:")
    print(wilcox.test(NDUFB7~Condition, data=df_cm))
    print(wilcox.test(Ferro~Condition, data=df_cm))
    print(wilcox.test(OXPHOS~Condition, data=df_cm))
  }
}

# 保存统计表
out <- "~/Projects/NDUFB7_HF_2026_04_20/03_results/02_tables/V68_GSE183852_Condition_Stats.csv"
write.csv(df, out, row.names=FALSE)
message("✅ 保存: ", out)
