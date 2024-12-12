## ANCOM-BC2 差异丰度分析脚本

# 安装和加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")
install.packages("readr")
install.packages("writexl")
install.packages("tidyverse")

library(ANCOMBC)
library(phyloseq)
library(readr)
library(writexl)
library(tidyverse)

# 加载 ANCOMBC2 函数文件
source("ANCOMBC2.R")

# 读取数据（注意使用分号作为分隔符）
otu_data <- read_delim("OTU_Table.csv", delim = ";", col_names = TRUE)
metadata <- read_delim("Metadata.csv", delim = ";", col_names = TRUE)

# 将 OTU 数据从 tibble 转换为 data.frame
otu_data <- as.data.frame(otu_data)

# 假设第一列是 OTU ID，其余列是样本的丰度值
rownames(otu_data) <- otu_data[[1]]   # 将 OTU ID 设为行名
otu_data <- otu_data[, -1]            # 删除 OTU ID 列，保留样本丰度信息

# 转换 OTU 矩阵并设置样本名称
otu_matrix <- as.matrix(otu_data)

# 确保 OTU 矩阵的列名为样本 ID（这些是样本名称）
colnames(otu_matrix) <- as.character(trimws(colnames(otu_matrix)))

# 确保元数据中 SampleID 列存在且与 OTU 矩阵一致
metadata <- as.data.frame(metadata)
rownames(metadata) <- as.character(trimws(metadata$SampleID))

# 获取共同样本 ID 并重新排序
common_samples <- intersect(colnames(otu_matrix), rownames(metadata))
otu_matrix <- otu_matrix[, common_samples]
metadata <- metadata[common_samples, ]

# 创建 phyloseq 对象
otu_table <- otu_table(otu_matrix, taxa_are_rows = TRUE)  # 注意 taxa_are_rows = TRUE，因为 OTU 是行
sample_data <- sample_data(metadata)
physeq <- phyloseq(otu_table, sample_data)

# 用户参数选择
cat("请选择分析类型: 1) 两组比较, 2) 多组比较\n")
analysis_type <- as.integer(readline(prompt = "输入 1 或 2: "))

cat("请选择模型类型: 1) 无协变量, 2) 包含协变量, 3) 包含随机因素\n")
model_type <- as.integer(readline(prompt = "输入 1, 2 或 3: "))

# 手动选择 p 值校正方法
cat("请选择 p 值校正方法: 1) holm, 2) bonferroni, 3) BH (FDR), 4) none\n")
p_adj_selection <- as.integer(readline(prompt = "输入 1, 2, 3 或 4: "))

p_adj_method <- switch(p_adj_selection,
                       "holm",
                       "bonferroni",
                       "BH",
                       "none")

# 手动设置 prv_cut, lib_cut 和 pseudo_sens 参数
cat("请输入 prv_cut 值 (0-1 之间, 推荐 0.10):\n")
prv_cut <- as.numeric(readline(prompt = "prv_cut: "))

cat("请输入 lib_cut 值 (推荐 1000):\n")
lib_cut <- as.numeric(readline(prompt = "lib_cut: "))

cat("是否开启 pseudo_sens (假数效应分析)? 1) 是, 2) 否\n")
pseudo_sens_selection <- as.integer(readline(prompt = "输入 1 或 2: "))

pseudo_sens <- ifelse(pseudo_sens_selection == 1, TRUE, FALSE)

cat("请选择线程数量 (推荐 1-4 之间):\n")
n_cl <- as.integer(readline(prompt = "线程数: "))

# 用户设置 iter_control 参数
cat("请输入迭代控制参数中的最大迭代次数 (max_iter),100: \n")
iter_max_iter <- as.integer(readline(prompt = "输入最大迭代次数: "))

cat("请输入迭代控制参数中的收敛公差 (tol),1e-5: \n")
iter_tol <- as.numeric(readline(prompt = "输入收敛公差: "))

# 构建公式
group_var <- "Group"  # 假设分组变量名为 'Group'
formula <- group_var

if (model_type == 2) {
  cat("请输入协变量名称，用空格分隔: ")
  covariates <- readline()
  formula <- paste(formula, covariates, sep = " + ")
} else if (model_type == 3) {
  cat("请输入随机因素名称: ")
  random_effect <- readline()
  formula <- paste(formula, "(1 |", random_effect, ")")
}

# 运行 ANCOM-BC2 分析
if (analysis_type == 1) {
  # 两组比较
  result <- ancombc2(
    data = physeq,            # 传递 phyloseq 对象
    fix_formula = formula,    # 线性模型公式
    group = group_var,        # 分组变量
    p_adj_method = p_adj_method,  # 使用用户选择的校正方法
    #prv_cut = prv_cut,        # 源象隐形切割
    #lib_cut = lib_cut,        # 库大小切割
    #pseudo_sens = pseudo_sens,# 假数效应分析
    #n_cl = n_cl,              # 线程数量
    #iter_control = list(tol = iter_tol, max_iter = iter_max_iter, verbose = FALSE)  # 用户设置的迭代控制参数
  )
} else {
  # 多组比较
  result <- ancombc2(
    data = physeq,            # 传递 phyloseq 对象
    fix_formula = formula,    # 线性模型公式
    group = group_var,        # 分组变量
    p_adj_method = p_adj_method,  # 使用用户选择的校正方法
    #prv_cut = prv_cut,        # 源象隐形切割
    #lib_cut = lib_cut,        # 库大小切割
    #pseudo_sens = pseudo_sens,# 假数效应分析
    #pairwise = TRUE,          # 启用两两比较
    #n_cl = n_cl,              # 线程数量
    #iter_control = list(tol = iter_tol, max_iter = iter_max_iter)  # 用户设置的迭代控制参数
  )
}

# 提取结果
res_prim <- result$res
res_pair <- result$res_pair

# 转化为数据表格格式
if (!is.data.frame(res_prim)) {
  res_prim <- as.data.frame(res_prim)
}

if (!is.data.frame(res_pair)) {
  res_pair <- as.data.frame(res_pair)
}

# 输出结果为 Excel 文件
write_xlsx(list(Primary_Analysis = res_prim, Pairwise_Comparisons = res_pair), "ancombc2_results.xlsx")

# 打印结果
print("分析完成。结果已保存为 ancombc2_results.xlsx")

# OTU_Table 示例格式
# OTU_Table 数据应该是一个包含 OTU ID 和样本丰度的矩阵，第一列是 OTU ID，第一行是样本名，剩余部分是丰度值。
# 示例：
# | OTU ID  | Sample1 | Sample2 | Sample3 |
# |---------|---------|---------|---------|
# | OTU1    |    5    |    8    |    0    |
# | OTU2    |    3    |    1    |    4    |
# | OTU3    |    0    |    2    |    7    |

# Metadata 示例格式
# Metadata 数据应包含每个样本的元数据信息，例如分组信息、协变量或随机因素。
# 示例：
# | SampleID | Group  | Age | Gender | Batch |
# |----------|--------|-----|--------|-------|
# | Sample1  | Group1 |  25 |   M    |   1   |
# | Sample2  | Group2 |  30 |   F    |   1   |
# | Sample3  | Group1 |  28 |   M    |   2   |

# 在不同参数选择下：
# 1. **两组比较，无协变量或随机因素**：只需包含 'SampleID' 和 'Group' 列。
# 2. **多组比较，包含协变量**：除了 'SampleID' 和 'Group'，还需要包含协变量列，例如 'Age'、'Gender'。
# 3. **包含随机因素**：需要包含 'SampleID'、'Group' 以及随机因素列，例如 'Batch'，表示实验批次或其他随机效应。

# 结果说明
# ANCOM-BC2 结果数据包含以下列：
# - **taxon**: 微生物分类单元 (例如物种或属)。
# - **lfc_(Intercept)** 和 **lfc_Group**: 每个分类单元的对数折叠变化 (log fold change)。
#   - `(Intercept)`: 代表基线组的对数折叠变化，表示在所有条件相同的情况下，微生物的基础丰度，可以理解为是对照组丰度。
#   - `Group`: 代表不同实验组与基线组之间的差异，表示某组与基线组相比的对数折叠变化。例如，`lfc_GroupRHCA` 表示 `GroupRHCA` 与基线组之间的差异。
#   - 实际生物学意义：
#     - 截距 (`(Intercept)`) 表示基线条件下的微生物丰度。基线条件是指没有受到任何处理或默认实验条件下的状态，它代表了在这种默认条件下微生物的自然丰度水平。这可以作为比较的基准，用于评估不同处理组如何改变微生物丰度。
#     - `group` 列则反映不同实验处理组与基线之间的相对丰度变化。显著的对数折叠变化意味着该组中的微生物丰度显著高于或低于基线组。
# - **se_(Intercept)** 和 **se_Group**: 对应的标准误差，表示估计的精度。
# - **W_(Intercept)** 和 **W_Group**: 统计量，等于 lfc / se，用于判断丰度是否显著不同于零。
# - **p_(Intercept)** 和 **p_Group**: p 值，基于双侧 Z 检验，用于检验丰度是否显著不同于零。
# - **q_(Intercept)** 和 **q_Group**: 调整后的 p 值，使用指定的多重比较校正方法。用于控制假阳性率。
# - **diff_(Intercept)** 和 **diff_Group**: 布尔值，指示该分类单元在该组是否显著 (TRUE 表示显著)。
# - **passed_ss_(Intercept)** 和 **passed_ss_Group**: 布尔值，指示该分类单元在灵敏度分析中是否通过。
#
# 生物学意义：
# - **基线组的丰度是否显著不同于零**：
#   - 如果基线组的丰度显著不同于零，表示该微生物在基线条件下普遍存在且具有重要的生态角色，能够作为对比参考点。
#   - 如果不显著不同于零，则该微生物在基线条件下稀少，可能对生态系统的影响有限，因此对其他组的比较意义也较弱。
# - **实验组与基线组的比较**：
#   - 通过比较实验组与基线组的对数折叠变化，可以评估实验处理是否对微生物的丰度产生显著影响，从而揭示不同条件下微生物群落结构的变化。
#
# 使用这些结果，您可以确定哪些微生物在不同分组之间表现出显著的丰度差异。对于多组比较，`res_pair` 将包含每组之间的对比结果。
