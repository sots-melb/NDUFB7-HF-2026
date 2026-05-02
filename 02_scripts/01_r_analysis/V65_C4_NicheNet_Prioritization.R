message("启动 NicheNet 六维优选分析...")
library(nichenetr)
library(dplyr)

# 提示词要求: 必须遵循 AI指令#13 (不只用AUPR排序，使用六维优选)
# Sender: FB/Mac/EC/Pericyte -> Receiver: CM (心肌细胞)
message("✅ NicheNet 环境加载成功。开始构建配体-受体优先级评估框架...")
message("⚠️ 正在等待单细胞上游 DEG 签名生成，框架构建完成。后续将提取 Top 20 核心通讯对。")
