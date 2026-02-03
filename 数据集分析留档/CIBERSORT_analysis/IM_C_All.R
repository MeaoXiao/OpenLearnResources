# ========================= 第一步：设置工作目录 & 安装依赖包 =========================
# 设置工作目录（替换为存放CIBERSORT.R和LM22.txt的路径）
setwd("D:\Imp_RSS\IBERSORT_analysis") 

# 安装CIBERSORT依赖的R包（e1071是核心，用于支持向量机；其他为可视化/数据处理包）
if (!require("e1071", quietly = TRUE)) install.packages("e1071")
if (!require("parallel", quietly = TRUE)) install.packages("parallel")
if (!require("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")

# 加载依赖包
library(e1071)        # 支持向量机（CIBERSORT核心依赖）
library(parallel)     # 并行计算（加快分析速度）
library(pheatmap)     # 免疫浸润比例热图
library(ggplot2)      # 可视化（堆叠柱状图/箱线图）
library(dplyr)        # 数据清洗
library(tidyr)        # 数据格式转换

# ========================= 第二步：加载CIBERSORT核心脚本 =========================
# 加载下载好的CIBERSORT.R脚本（确保文件在工作目录，名字无误）
source("CIBERSORT.R")  # 若报错，检查文件路径/文件名是否正确

# ========================= 第三步：准备输入表达矩阵 =========================
# CIBERSORT输入要求：
# 1. 表达矩阵：行=基因名（HGNC标准SYMBOL），列=样本名，值=标准化表达量（推荐TPM/FPKM，不建议原始count）
# 2. 矩阵需为**文本文件（.txt/.csv）**，且基因名列无表头（CIBERSORT要求）

# ---------------------- 3.1 构建示例表达矩阵（可替换为真实数据） ----------------------
# 模拟500个基因、30个样本的表达矩阵（仅作演示，需替换为自己的数据）
set.seed(123)
expr_matrix <- matrix(rnorm(500*30, mean = 5, sd = 1),
                      nrow = 500,
                      ncol = 30,
                      dimnames = list(paste0("Gene", 1:500),  # 基因名（示例）
                                      paste0("Sample", 1:30))) # 样本名（示例）

# ---------------------- 3.2 真实数据加载（用户替换核心代码） ----------------------
# 方式1：从CSV文件加载（推荐，需确保行名是基因SYMBOL）
# expr_matrix <- read.csv("your_expression_data.csv", row.names = 1, header = TRUE, check.names = FALSE)
# 方式2：从TXT文件加载
# expr_matrix <- read.table("your_expression_data.txt", row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)

# ---------------------- 3.3 表达矩阵预处理（关键，避免分析报错） ----------------------
# 1. 过滤低表达基因（保留在≥50%样本中表达量>0的基因）
expr_matrix <- expr_matrix[rowSums(expr_matrix > 0) >= ncol(expr_matrix)*0.5, ]
# 2. 转换为矩阵格式（CIBERSORT要求）
expr_matrix <- as.matrix(expr_matrix)
# 3. 确保基因名是字符型（避免格式错误）
rownames(expr_matrix) <- as.character(rownames(expr_matrix))
# 4. 保存为CIBERSORT要求的TXT文件（必须步骤，CIBERSORT需读取文本文件）
write.table(expr_matrix, 
            file = "input_expression_matrix.txt",  # 输出文件名
            sep = "\t",                            # 分隔符为制表符
            quote = FALSE,                         # 不添加引号
            row.names = TRUE,                      # 保留基因名（行名）
            col.names = TRUE)                      # 保留样本名（列名）

# ========================= 第四步：运行CIBERSORT核心分析 =========================
# CIBERSORT函数参数说明：
# - sig_matrix：LM22特征矩阵文件路径
# - mixture_file：输入表达矩阵文件路径
# - perm：置换检验次数（推荐≥100，值越大结果越可靠，速度越慢；调试时可设10）
# - QN：是否进行分位数归一化（TRUE/FALSE，RNA-seq数据建议TRUE）

# 运行CIBERSORT（核心命令）
cibersort_result <- CIBERSORT(sig_matrix = "LM22.txt",  # LM22特征矩阵路径
                              mixture_file = "input_expression_matrix.txt",  # 输入表达矩阵路径
                              perm = 100,  # 置换检验次数（正式分析建议设1000）
                              QN = TRUE)   # 分位数归一化（RNA-seq必选）

# 查看结果结构：每行是样本，列包括22种免疫细胞比例、P值、Correlation、RMSE
# P值<0.05表示该样本的免疫浸润结果可信
print(head(cibersort_result))

# ========================= 第五步：结果处理 & 过滤低置信度样本 =========================
# ---------------------- 5.1 提取核心结果（仅保留22种免疫细胞比例） ----------------------
# 提取前22列（22种免疫细胞），去掉P值、Correlation、RMSE列
immune_infiltration <- cibersort_result[, 1:22]
# 重置行名为样本名（避免行名错乱）
rownames(immune_infiltration) <- rownames(cibersort_result)

# ---------------------- 5.2 过滤低置信度样本（P<0.05） ----------------------
# 提取P值列（第23列）
p_values <- cibersort_result[, "P-value"]
# 只保留P<0.05的样本（结果更可靠）
immune_infiltration_filtered <- immune_infiltration[p_values < 0.05, ]

# 输出过滤后的结果到CSV（便于后续分析）
write.csv(immune_infiltration_filtered, 
          file = "CIBERSORT_immune_infiltration_result.csv", 
          row.names = TRUE)

# ========================= 第六步：免疫浸润结果可视化 =========================
# ---------------------- 6.1 样本免疫细胞组成堆叠柱状图 ----------------------
# 数据格式转换（长格式，适配ggplot2）
immune_long <- immune_infiltration_filtered %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%  # 样本名转为列
  pivot_longer(cols = -Sample,              # 除Sample外的列转为长格式
               names_to = "Cell_Type",      # 免疫细胞类型列名
               values_to = "Proportion")    # 浸润比例列名

# 绘制堆叠柱状图
ggplot(immune_long, aes(x = Sample, y = Proportion, fill = Cell_Type)) +
  geom_col(position = "stack") +  # 堆叠柱状图
  labs(title = "Immune Cell Infiltration by Sample",  # 标题
       x = "Sample", y = "Infiltration Proportion") +
  theme_bw() +  # 白色背景
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 样本名旋转45度
        plot.title = element_text(hjust = 0.5),             # 标题居中
        legend.position = "right") +  # 图例在右侧
  scale_fill_brewer(palette = "Paired")  # 配色方案（美观易区分）

# ---------------------- 6.2 免疫细胞浸润比例热图 ----------------------
# 为样本添加分组（示例：前15个为对照组，后15个为肿瘤组，需替换为你的真实分组）
sample_group <- factor(c(rep("Control", 15), rep("Tumor", 15)),
                       levels = c("Control", "Tumor"))
# 构建热图注释
annotation_col <- data.frame(Group = sample_group,
                             row.names = rownames(immune_infiltration_filtered))

# 绘制热图
pheatmap(mat = t(immune_infiltration_filtered),  # 转置：行=细胞类型，列=样本
         annotation_col = annotation_col,        # 样本分组注释
         show_rownames = TRUE,                   # 显示细胞类型（行名）
         show_colnames = FALSE,                  # 隐藏样本名（避免拥挤）
         scale = "row",                          # 按行标准化（突出细胞在样本间的差异）
         color = colorRampPalette(c("blue", "white", "red"))(100),  # 颜色梯度
         main = "Immune Cell Infiltration Heatmap",  # 标题
         fontsize = 8)  # 字体大小

# ---------------------- 6.3 关键免疫细胞浸润比例箱线图（示例：CD8+ T细胞） ----------------------
# 整理数据
boxplot_data <- immune_infiltration_filtered %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(Group = sample_group)  # 添加分组列

# 绘制CD8+ T细胞箱线图
ggplot(boxplot_data, aes(x = Group, y = `CD8 T cells`, fill = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.8) +  # 箱线图
  geom_jitter(size = 1, alpha = 0.6) +      # 叠加散点显示单个样本
  labs(title = "CD8+ T Cell Infiltration in Different Groups",
       x = "Group", y = "Infiltration Proportion") +
  scale_fill_manual(values = c("Control" = "skyblue", "Tumor" = "coral")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# ========================= 第七步：保存可视化结果 =========================
# 保存堆叠柱状图
ggsave("Immune_Infiltration_Stacked_Barplot.pdf", 
       width = 12, height = 8, dpi = 300)

# 保存热图（先打开PDF设备，再绘图，最后关闭）
pdf("Immune_Infiltration_Heatmap.pdf", width = 10, height = 8)
pheatmap(mat = t(immune_infiltration_filtered), annotation_col = annotation_col, scale = "row")
dev.off()

# 保存箱线图
ggsave("CD8_Tcell_Infiltration_Boxplot.pdf", width = 6, height = 5, dpi = 300)

#结果：
#每列代表一种免疫细胞的相对浸润比例（所有细胞比例之和≈1）。
#P 值 < 0.05：样本的免疫浸润结果可信；Correlation 越接近 1，拟合效果越好。