# ========================= 第一步：安装必要的R包 =========================
# GSVA属于Bioconductor包，需用BiocManager安装（先安装BiocManager）
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")  # 安装BiocManager（CRAN包）
}

# 安装核心包和辅助包（按需安装，已安装则跳过）
BiocManager::install(c("GSVA",        # GSVA核心分析包
                       "msigdbr",     # 获取MSigDB基因集（替代clusterProfiler的getGSEA函数）
                       "org.Hs.eg.db",# 人类基因注释包（用于基因名转换）
                       "pheatmap",    # 热图可视化
                       "ggplot2",     # 通用可视化
                       "dplyr"),      # 数据处理
                     update = FALSE,  # 不更新已有包（加快速度）
                     ask = FALSE)     # 无需手动确认

# ========================= 第二步：加载安装好的R包 =========================
library(GSVA)          # 加载GSVA核心包
library(msigdbr)       # 加载基因集获取包
library(org.Hs.eg.db)  # 加载人类基因注释数据库
library(pheatmap)      # 加载热图绘制包
library(ggplot2)       # 加载绘图包
library(dplyr)         # 加载数据处理包

# ========================= 第三步：准备输入数据 =========================
# GSVA需要两个核心输入：
# 1. 基因表达矩阵（行=基因名，列=样本名，值=表达量，推荐log2标准化后的数据）
# 2. 基因集（格式：列表，每个元素是一个基因集，元素值为该基因集的基因名）

# ---------------------- 3.1 构建示例表达矩阵（可替换为自己的真实数据） ----------------------
# 模拟1000个基因、20个样本的表达矩阵（用户需替换为自己的表达矩阵，如从CSV/TCGA下载）
set.seed(123)  # 设置随机种子，保证结果可重复
expr_matrix <- matrix(rnorm(1000*20),  # 正态分布随机数
                      nrow = 1000,     # 1000个基因
                      ncol = 20,       # 20个样本
                      dimnames = list(paste0("Gene", 1:1000),  # 基因名（示例）
                                      paste0("Sample", 1:20))) # 样本名（示例）

# 真实数据加载示例（用户替换）：
# expr_matrix <- read.csv("your_expression_matrix.csv", row.names = 1, header = TRUE)
# 注意：表达矩阵需满足：行是基因名（SYMBOL/ENTREZ ID），列是样本，值为标准化表达量

# ---------------------- 3.2 过滤低表达基因（可选但推荐） ----------------------
# 过滤掉在超过50%样本中表达量为0/极低的基因（减少噪音）
expr_matrix <- expr_matrix[rowSums(expr_matrix > 0) > ncol(expr_matrix)*0.5, ]

# ---------------------- 3.3 获取基因集（以MSigDB的H hallmark基因集为例） ----------------------
# 使用msigdbr获取人类H hallmark基因集（最常用的功能基因集）
msig_df <- msigdbr(species = "Homo sapiens",  # 物种：人类（小鼠用"Mus musculus"）
                   category = "H")            # 基因集类别：H（hallmark），也可选C1-C7

# 将基因集转换为GSVA要求的列表格式（列表名=基因集名称，列表值=基因名）
gene_sets <- split(msig_df$gene_symbol,  # 基因名（SYMBOL）
                   msig_df$gs_name)     # 基因集名称

# 可选：筛选基因集（只保留在表达矩阵中存在的基因，提升分析效率）
gene_sets <- lapply(gene_sets, function(genes) {
  intersect(genes, rownames(expr_matrix))  # 只保留表达矩阵中存在的基因
})
# 移除空基因集（避免分析报错）
gene_sets <- gene_sets[sapply(gene_sets, length) > 0]

# ========================= 第四步：运行GSVA分析 =========================
# gsva()函数核心参数说明：
# - expr: 表达矩阵（行=基因，列=样本）
# - gset.idx.list: 基因集列表
# - method: GSVA计算方法（默认"gsva"，可选"ssgsea"/"zscore"/"plage"）
# - kcdf: 表达数据类型（"Gaussian"适用于微阵列/标准化RNA-seq；"Poisson"适用于原始count）
# - min.sz/max.sz: 基因集最小/最大基因数（过滤极端基因集）

gsva_result <- gsva(expr = expr_matrix,
                    gset.idx.list = gene_sets,
                    method = "gsva",        # 使用标准GSVA算法
                    kcdf = "Gaussian",      # 模拟数据为正态分布，真实RNA-seq标准化后也用这个
                    min.sz = 5,             # 基因集至少包含5个基因
                    max.sz = 500,           # 基因集最多包含500个基因
                    verbose = TRUE)         # 显示运行过程

# 查看GSVA结果：行=基因集，列=样本，值=GSVA得分（反映该基因集在样本中的激活程度）
print(head(gsva_result))

# ========================= 第五步：GSVA结果可视化 =========================
# ---------------------- 5.1 绘制GSVA得分热图（核心可视化） ----------------------
# 为样本添加分组（示例：前10个为对照组，后10个为处理组）
sample_group <- factor(c(rep("Control", 10), rep("Treatment", 10)),
                       levels = c("Control", "Treatment"))
# 构建热图注释（样本分组）
annotation_col <- data.frame(Group = sample_group,
                             row.names = colnames(gsva_result))

# 绘制热图
pheatmap(mat = gsva_result,                # GSVA得分矩阵
         annotation_col = annotation_col,  # 样本分组注释
         show_rownames = TRUE,             # 显示基因集名称（行名）
         show_colnames = FALSE,            # 隐藏样本名（避免拥挤）
         scale = "row",                    # 按行标准化（突出基因集在不同样本的差异）
         color = colorRampPalette(c("blue", "white", "red"))(100),  # 颜色梯度
         main = "GSVA Score Heatmap",      # 标题
         fontsize = 8)                     # 字体大小

# ---------------------- 5.2 绘制单个基因集得分箱线图（示例） ----------------------
# 选择第一个基因集（HALLMARK_TNFA_SIGNALING_VIA_NFKB）绘制箱线图
gs_name <- names(gene_sets)[1]
# 整理数据
boxplot_data <- data.frame(
  GSVA_Score = gsva_result[gs_name, ],    # 该基因集的GSVA得分
  Group = sample_group                    # 样本分组
)

# 绘制箱线图
ggplot(boxplot_data, aes(x = Group, y = GSVA_Score, fill = Group)) +
  geom_boxplot(width = 0.5, alpha = 0.8) +  # 箱线图
  geom_jitter(size = 1, alpha = 0.6) +      # 散点叠加
  labs(title = paste("GSVA Score of", gs_name),  # 标题
       x = "Group", y = "GSVA Score") +
  scale_fill_manual(values = c("Control" = "skyblue", "Treatment" = "coral")) +
  theme_bw() +  # 白色背景
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中

# ========================= 第六步：保存结果 =========================
# 保存GSVA得分矩阵到CSV文件
write.csv(gsva_result, "GSVA_Result.csv", row.names = TRUE)

# 保存热图和箱线图（PDF格式，高清）
pdf("GSVA_Heatmap.pdf", width = 10, height = 8)
pheatmap(mat = gsva_result, annotation_col = annotation_col, show_rownames = TRUE, show_colnames = FALSE, scale = "row")
dev.off()

pdf("GSVA_Boxplot.pdf", width = 6, height = 5)
print(ggplot(boxplot_data, aes(x = Group, y = GSVA_Score, fill = Group)) + geom_boxplot() + geom_jitter() + theme_bw())
dev.off()

#结果：
#GSVA 得分越高，代表该基因集在对应样本中越活跃；
#热图可直观展示基因集在不同分组样本的激活模式，
#箱线图可比较单个基因集在分组间的差异。