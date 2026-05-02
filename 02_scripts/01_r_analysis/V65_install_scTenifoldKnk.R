# R Script for scTenifoldKnk safe installation
message("开始安装 scTenifoldKnk 及其依赖...")

if(!require("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org/")
}

if(!require("scTenifoldKnk", quietly = TRUE)) {
    message("尝试从 GitHub 安装 scTenifoldKnk...")
    tryCatch({
        devtools::install_github("SystemsBiologylmu/scTenifoldKnk")
        message("✅ scTenifoldKnk 安装成功!")
    }, error = function(e) {
        message("❌ 安装失败，错误信息: ", e$message)
        message("替代方案: 如果私有库失效，请在后续分析中降级使用 scTenifoldNet。")
    })
} else {
    message("✅ scTenifoldKnk 已安装，跳过。")
}
