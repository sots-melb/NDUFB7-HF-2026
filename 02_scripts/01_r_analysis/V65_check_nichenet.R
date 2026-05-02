if(!require("nichenetr", quietly = TRUE)) {
    message("正在安装 nichenetr...")
    devtools::install_github("saeyslab/nichenetr")
} else {
    message("✅ nichenetr 已安装。")
}
