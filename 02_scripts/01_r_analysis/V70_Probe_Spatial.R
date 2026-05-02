suppressMessages(library(GEOquery))

gse_id <- "GSE290094"
message("▶ 正在从 GEO 获取 ", gse_id, " 的元数据 (您上传附件中的小鼠 Visium 数据)...")

tryCatch({
  gse <- getGEO(gse_id, GSEMatrix=TRUE, getGPL=FALSE, AnnotGPL=FALSE)
  if(length(gse) > 0) {
    pdata <- pData(gse[[1]])
    message("✅ 数据集访问成功！")
    message("   样本总数: ", nrow(pdata))
    message("   包含的组织类型/条件:")
    
    # 打印前几个样本的核心特征，查看是否有假手术组 (Sham) 和缺血组
    char_cols <- grep("characteristics_ch1", colnames(pdata), value=TRUE)
    if(length(char_cols) > 0) {
      for(i in 1:min(5, nrow(pdata))) {
        message("   Sample ", i, ": ", paste(pdata[i, char_cols], collapse=" | "))
      }
    }
  }
}, error = function(e) {
  message("❌ 获取元数据失败，可能网络连接超时: ", conditionMessage(e))
})
