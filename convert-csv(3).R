# 安装并加载必要的包
if (!require("readxl")) install.packages("readxl")
if (!require("writexl")) install.packages("writexl")

library(readxl)

# 读取Excel文件并将其转换为分号分隔符的CSV文件
convert_excel_to_csv <- function(input_file, output_file, sheet = 1) {
  # 读取Excel文件
  data <- read_excel(input_file, sheet = sheet)
  
  # 写入CSV文件，分号分隔符
  write.table(data, file = output_file, sep = ";", row.names = FALSE, col.names = TRUE, quote = TRUE)
}
# 使用示例
# 将"your_excel_file.xlsx"替换为你的Excel文件名
# 将"output_file.csv"替换为你想要输出的CSV文件名
convert_excel_to_csv("Metadata.xlsx", "Metadata.csv")
