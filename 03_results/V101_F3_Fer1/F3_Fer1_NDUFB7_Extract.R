
# F3_Fer1_NDUFB7_Extract.R
# 从GSE243655 series matrix提取NDUFB7并做配对检验

library(GEOquery)
gse <- getGEO("GSE243655", GSEMatrix = TRUE)
exprs <- exprs(gse[[1]])  # 表达矩阵

# 找NDUFB7探针
probes <- rownames(exprs)
ndufb7_probe <- probes[grep("NDUFB7", probes, ignore.case = TRUE)][1]

if (!is.na(ndufb7_probe)) {
  ndufb7_vals <- exprs[ndufb7_probe, ]
  
  # 读取分组表
  grp <- read.csv("V101_GSE243655_groups_confirmed.csv")
  
  # 配对检验
  dmso <- ndufb7_vals[grp$Treatment == "DMSO"]
  fer1 <- ndufb7_vals[grp$Treatment == "Fer-1"]
  
  # 配对t-test
  tt <- t.test(dmso, fer1, paired = TRUE)
  
  # 输出
  cat("NDUFB7 Probe:", ndufb7_probe, "\n")
  cat("DMSO mean:", mean(dmso), "\n")
  cat("Fer-1 mean:", mean(fer1), "\n")
  cat("Delta:", mean(fer1 - dmso), "\n")
  cat("Paired t-test p:", tt$p.value, "\n")
  
  # 保存
  res <- data.frame(
    Patient = 1:4,
    DMSO = dmso,
    Fer1 = fer1,
    Delta = fer1 - dmso
  )
  write.csv(res, "V101_Fer1_NDUFB7_paired.csv", row.names = FALSE)
}

