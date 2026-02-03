单细胞分析基本流程（Seurat版本5.0.0，以GSE188545为例）
# ==================== [T] 0.加载包部分 ====================
##1. 环境准备与数据加载
# 先安装所有缺失的包和依赖
# 安装clustree包（从CRAN下载）
if (!require("clustree", quietly = TRUE)) {
  install.packages("clustree", dependencies = TRUE)
}

# 安装ROCR包（Seurat的依赖包）
if (!require("ROCR", quietly = TRUE)) {
  install.packages("ROCR", dependencies = TRUE)
}

# 安装Seurat包（如果之前安装不完整，重新安装确保依赖齐全）
if (!require("Seurat", quietly = TRUE)) {
  install.packages("Seurat", dependencies = TRUE)
}

# 安装其他已加载成功但确保版本兼容的包（可选，确保环境一致）
required_packages <- c("ggplot2", "dplyr", "gridExtra")
for (pkg in required_packages) {
  if (!require(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

##1.1 加载必要的R包

# 加载各种分析所需的R包
library(clustree)    # 用于聚类树可视化，帮助选择最佳分辨率
library(Seurat)      # 单细胞分析的核心包，提供主要分析功能
library(ggplot2)     # 数据可视化包，用于创建高质量的图形
library(dplyr)       # 数据操作包，提供数据筛选、汇总等功能
library(gridExtra)   # 图形排列包，用于组合多个图形
library(Seurat)      # 再次加载确保环境正确
Read10X()
#如果报没有该函数则Seurat加载失败
# ==================== [T] 1. 工作路径以及文件预处理 ====================
#1.2 设置工作目录并读取数据

# 定义存储原始数据的目录路径
dir <- "D:/Imp_RSS/GSE188545_RAW"

# 获取目录中所有样本文件的名称
samples <- list.files(dir, pattern = "\\.gz$")
# 验证
print(samples)
print(length(samples))

# 提取每个文件名的前10个字符作为样本简称
samples_short <- substr(samples, 1, 10)
unique_samples_short <- unique(samples_short)
# 获取唯一样本ID
unique_samples <- unique(samples_short)
cat("发现", length(unique_samples), "个唯一样本:\n")
print(unique_samples)

# 批量读取所有样本数据并创建Seurat对象列表
sceList <- lapply(unique_samples, function(sample_id){
  cat("\n=== 处理样本:", sample_id, "===\n")
  
  # 1. 找到该样本的三个文件
  sample_files <- samples[grepl(paste0("^", sample_id, "_"), samples)]
  
  barcode_file <- sample_files[grepl("barcodes", sample_files)]
  gene_file <- sample_files[grepl("genes", sample_files)]
  matrix_file <- sample_files[grepl("matrix", sample_files)]
  
  cat("找到文件:\n")
  cat("  barcodes:", barcode_file, "\n")
  cat("  genes:", gene_file, "\n")
  cat("  matrix:", matrix_file, "\n")
  
  # 2. 创建临时目录
  temp_dir <- tempfile(pattern = paste0(sample_id, "_"))
  dir.create(temp_dir)
  cat("临时目录:", temp_dir, "\n")
  
  # 3. 复制并重命名文件为标准名称
  file.copy(file.path(dir, barcode_file), 
            file.path(temp_dir, "barcodes.tsv.gz"))
  file.copy(file.path(dir, gene_file), 
            file.path(temp_dir, "features.tsv.gz"))
  file.copy(file.path(dir, matrix_file), 
            file.path(temp_dir, "matrix.mtx.gz"))
  
  # 4. 使用Read10X读取数据
  tryCatch({
    counts <- Seurat::Read10X(data.dir = temp_dir)
    
    # 5. 创建Seurat对象
    seurat_obj <- Seurat::CreateSeuratObject(
      counts = counts,
      project = sample_id,
      min.cells = 5,
      min.features = 200
    )
    
    # 6. 清理临时目录
    unlink(temp_dir, recursive = TRUE)
    
    cat("✓ 成功创建: ", ncol(seurat_obj), "个细胞\n")
    
    return(seurat_obj)
    
  }, error = function(e) {
    # 清理临时目录
    unlink(temp_dir, recursive = TRUE)
    cat("✗ 读取失败:", e$message, "\n")
    return(NULL)
  })
})
# 移除失败的项目
sceList <- sceList[!sapply(sceList, is.null)]

# 命名列表
names(sceList) <- unique_samples[!sapply(sceList, is.null)]

# 总结结果
cat("成功处理:", length(sceList), "/", length(unique_samples), "个样本\n")

if (length(sceList) > 0) {
  cat("成功样本:", names(sceList), "\n")
  
  # 显示第一个样本的信息
  cat("\n第一个样本详情:\n")
  print(sceList[[1]])
}

#1.3 合并所有样本数据

# 查看样本简称
samples_short

# 合并所有Seurat对象，为每个样本的细胞添加前缀标识
merged_seurat <- merge(x = sceList[[1]],
                       y = sceList[-1],
                       add.cell.ids = unique_samples_short)
# 查看合并后的元数据前几行
head(merged_seurat@meta.data)

# 查看数据维度（基因数×细胞数）
dim(merged_seurat)

# 统计各样本的细胞数量
table(merged_seurat$orig.ident)
#1.4 添加样本信息

# 查看当前orig.ident的值
cat("当前的orig.ident值示例:\n")
print(head(merged_seurat$orig.ident, 10))

# 问题：orig.ident只是GSM编号，不包含AD/HC信息
# 所以需要从原始文件名中提取AD/HC信息
# 可以先创建一个样本信息映射表
sample_info <- data.frame(
  sample_id = c("GSM5685287", "GSM5685288", "GSM5685289", 
                "GSM5685290", "GSM5685291", "GSM5685292",
                "GSM5685293", "GSM5685294", "GSM5685295",
                "GSM5685296", "GSM5685297", "GSM5685298"),
  condition = c("AD", "AD", "AD", "AD", "AD", "AD",
                "HC", "HC", "HC", "HC", "HC", "HC"),
  sample_type = c("AD", "AD", "AD", "AD", "AD", "AD",
                  "HC", "HC", "HC", "HC", "HC", "HC")
)

# 载入映射
merged_seurat$condition <- ifelse(
  merged_seurat$orig.ident %in% c("GSM5685287", "GSM5685288", "GSM5685289",
                                  "GSM5685290", "GSM5685291", "GSM5685292"),
  "AD",
  "HC"
)
merged_seurat$sample_type <- merged_seurat$condition
# 检验
cat("\n=== 疾病状态分布 ===\n")
table(merged_seurat$condition)
print(table(merged_seurat$condition))

cat("\nAD样本数:", sum(merged_seurat$condition == "AD"), "\n")
cat("HC样本数:", sum(merged_seurat$condition == "HC"), "\n")

merged_seurat$sample <- merged_seurat$orig.ident
# 查看疾病状态分布
table(merged_seurat$condition)
# ==================== [T] 2. 质控指标计算与可视化 ====================
#2. 数据质控
#2.1 计算质控指标

# 计算线粒体基因比例（细胞应激指标）
merged_seurat[["percent_mito"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

# 计算血红蛋白基因比例（红细胞污染指标）
merged_seurat[["percent_hb"]] <- PercentageFeatureSet(merged_seurat, pattern = "^HB[ABDGQEMZ]")

# 计算核糖体基因比例（蛋白质合成活性指标）
merged_seurat[["percent_ribo"]] <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]")

#2.2 质控指标可视化

# 先检查数据
cat("检查可用metadata列:\n")
print(colnames(merged_seurat@meta.data))
library(ggplot2)
library(patchwork)
# 绘制质控指标的小提琴图
p1_qc <- VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0)
ggsave("545/1.QC_Violin.pdf", plot = p1_qc, width = 12, height = 5)
print(p1_qc)

# 创建质控指标散点图
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent_mito")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2_qc_scatter <- plot1 + plot2
ggsave("545/2.QC_Scatter.pdf", plot = p2_qc_scatter, width = 10, height = 5)
print(p2_qc_scatter)


# 图片解释：
# •	1.QC_Violin.pdf: 显示每个细胞的基因数(nFeature_RNA)、UMI数(nCount_RNA)和线粒体基因百分比分布。用于识别低质量细胞。
# •	2.QC_Scatter.pdf: 左图显示UMI数与线粒体基因百分比的关系，右图显示UMI数与基因数的关系，用于检测异常细胞。

#2.3 样本间质控比较

# 按样本分组绘制质量指标的小提琴图
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p4_sample_metrics <- VlnPlot(merged_seurat, group.by = "sample", features = feats, pt.size = 0, ncol = 3, same.y.lims = T) +
  scale_y_continuous(breaks = seq(0, 100, 5)) + NoLegend()
w <- length(unique(merged_seurat$orig.ident))/2 + 5
ggsave("545/4.Sample_QC_Metrics.pdf", plot = p4_sample_metrics, width = w, height = 5)
print(p4_sample_metrics)
# 图片解释：
# •	4.Sample_QC_Metrics.pdf: 比较不同样本间的线粒体基因、核糖体基因和血红蛋白基因百分比，识别样本特异性质量问题。

# 2.4 细胞过滤

# 基于质控指标过滤低质量细胞
merged_seurat<- JoinLayers(merged_seurat)
merged_seurat <- subset(
  merged_seurat,
  subset = nFeature_RNA > 200 &     # 保留基因数在200-5000之间的细胞
    nFeature_RNA < 5000 &
    percent_mito < 5                  # 保留线粒体基因百分比小于5%的细胞
)
#如以下代码报错，先运行merged_seurat<- JoinLayers(merged_seurat)


# ==================== [T] 3. 双细胞检测与去除 ====================
# ==================== [T] 3.1 方法一：DoubletFinder ====================
# ==================== [T] 3.1.1 DoubletFinder准备步骤 ====================
# 从GitHub安装DoubletFinder包
devtools::install_github("chris-mcginnis-ucsf/DoubletFinder", force = TRUE)

# 设置预期的双细胞率
doublet_rate <- 0.016  # 1.6%的双细胞率
library(DoubletFinder)
seurat.data <- merged_seurat
print(seurat.data)
# 标准预处理流程
seurat.data <- NormalizeData(seurat.data)
seurat.data <- FindVariableFeatures(seurat.data, selection.method = "vst", nfeatures = 2000)
seurat.data <- ScaleData(seurat.data)
# 这一步不降维FindNeighbors无法正确运行
seurat.data <- RunPCA(seurat.data)
# 检查PCA结果，确定使用的主成分数量
# 可视化肘部图，帮助判断主成分“拐点”
ElbowPlot(seurat.data, ndims = 50)
# 通常会选择肘部图拐点之后的PCs，例如1:30或1:40
seurat.data <- FindNeighbors(seurat.data, reduction = "pca", dims = 1:50)
seurat.data <- FindClusters(seurat.data, resolution = 0.7)

# 降维可视化
seurat.data <- seurat.data %>%
  RunPCA(npcs = 50, verbose = F) %>%
  RunUMAP(reduction = "pca", dims = 1:50, verbose = F)

# 参数扫描寻找最佳pK值
sweep.res.list <- paramSweep(seurat.data, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
# ==================== [F] 3.1.2 DoubletFinder直接处理方案 ====================
# 计算预期的双细胞数并运行DoubletFinder
nExp_poi <- round(doublet_rate * ncol(seurat.data))
seurat.data <- doubletFinder(seurat.data, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, sct = FALSE)

# 过滤掉双细胞，只保留单细胞
seurat.data <- subset(seurat.data, DF.classifications_0.25_0.005_105 == "Singlet")
# ==================== [T] 3.1.3 DoubletFinder循环处理方案 ====================
# 循环处理法
if (!"sample" %in% colnames(seurat.data@meta.data)) {
  seurat.data$sample <- seurat.data$orig.ident
}

# 获取所有样本名称
samples <- unique(seurat.data$sample)
doublet_results <- list() # 用于存储每个样本的双细胞预测列

# 对每个样本循环处理
for (samp in samples) {
  cat("\n正在处理样本:", samp, "\n")
  
  # 1. 提取单个样本的子集
  sub_seurat <- subset(seurat.data, subset = sample == samp)
  cat("  细胞数:", ncol(sub_seurat), "\n")
  
  # 2. 对该样本进行标准预处理
  sub_seurat <- NormalizeData(sub_seurat)
  sub_seurat <- FindVariableFeatures(sub_seurat, selection.method = "vst", nfeatures = 2000)
  sub_seurat <- ScaleData(sub_seurat)
  sub_seurat <- RunPCA(sub_seurat)
  
  # 3. 为DoubletFinder准备（聚类）
  sub_seurat <- FindNeighbors(sub_seurat, dims = 1:20)
  sub_seurat <- FindClusters(sub_seurat, resolution = 0.7)
  
  # 4. 优化参数并运行DoubletFinder（样本单独计算）
  sweep.res <- paramSweep(sub_seurat, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  nCells_sample <- ncol(sub_seurat)
  nExp_sample <- round(nCells_sample * nCells_sample / 100000 * 0.8) # 每1000细胞0.8%的近似计算
  
  # 5. 运行DoubletFinder（内存需求大大降低）
  sub_seurat <- doubletFinder(sub_seurat,
                                 PCs = 1:20,
                                 pN = 0.25,
                                 pK = optimal_pK,
                                 nExp = nExp_sample,
                                 sct = FALSE)
  
  # 6. 提取预测结果
  df_column_name <- grep("^DF.classifications", colnames(sub_seurat@meta.data), value = TRUE)[1]
  doublet_results[[samp]] <- data.frame(
    Cell = colnames(sub_seurat),
    Doublet_Label = sub_seurat@meta.data[[df_column_name]],
    Sample = samp
  )
}

# 合并所有样本的结果
library(dplyr)
doublet_df <- bind_rows(doublet_results)

# 将结果添加到原始seurat.data对象中
seurat.data$Doublet_Label <- "Singlet" # 默认值
rownames(doublet_df) <- doublet_df$Cell
matching_cells <- intersect(rownames(doublet_df), colnames(seurat.data))
seurat.data$Doublet_Label[matching_cells] <- doublet_df[matching_cells, "Doublet_Label"]

# 查看总体结果
cat("\n=== 双细胞预测汇总 ===\n")
table(seurat.data$Doublet_Label)
table(seurat.data$sample, seurat.data$Doublet_Label)
seurat.data <- RunTSNE(seurat.data, reduction = "pca", dims = 1:20)
p <- DimPlot(seurat.data,
             reduction = "tsne",
             group.by = "Doublet_Label",       # 使用存储预测结果的列
             cols = c("Singlet" = "lightgrey", "Doublet" = "red"), # 自定义颜色
             pt.size = 0.5,                    # 点的大小，可根据细胞数量调整
             order = c("Doublet"))              # 将双细胞点绘制在最上层，避免被遮盖
print(p)

#图片解释：
#•	5.DoubletFinder_Detection.pdf: 显示DoubletFinder识别的双细胞（doublets）和单细胞（singlets）在t-SNE图中的分布。
# ==================== [F] 3.2 方法二：scDblFinder ====================
#r
# 设置随机种子保证结果可重复
the.seed <- 1337L

# 使用scDblFinder检测双细胞
BiocManager::install("scDblFinder",force = TRUE)
library("scDblFinder")
merged_seurat$scDblFinder.class <- Seurat::as.Seurat(
  scDblFinder::scDblFinder(Seurat::as.SingleCellExperiment(merged_seurat))
)$scDblFinder.class

merged_seurat$scDblFinder.class <- unname(merged_seurat$scDblFinder.class == "doublet")

# 输出scDblFinder检测结果
print('scDblFinder doublets :')
print(table(merged_seurat$scDblFinder.class))
# ==================== [F] 3.3 方法三：scds混合方法 ====================
#r
BiocManager::install("scds")
library(scds)
# 使用scds包的混合方法检测双细胞
merged_seurat$scds.score <- scds::cxds_bcds_hybrid(
  Seurat::as.SingleCellExperiment(merged_seurat)
)$hybrid_score

merged_seurat$scds.class <- unname(merged_seurat$scds.score > 1)

# 输出scds检测结果
print('scds-hybrid doublets :')
print(table(merged_seurat$scds.class))
# ==================== [F] 3.4 双细胞检测结果分析与过滤原案 ====================

# 整合两种方法的双细胞检测结果（取并集）
merged_seurat$doublets_consensus.class <- merged_seurat$scDblFinder.class | merged_seurat$scds.class

# 输出共识双细胞结果
print('Consensus doublets :')
print(table(merged_seurat$doublets_consensus.class))

# 分析双细胞与细胞周期的关系
df_S.cycle_doublets <- data.frame(Phase = merged_seurat$Phase, doublets = merged_seurat$doublets_consensus.class)
print(table(df_S.cycle_doublets))

# 过滤掉双细胞
merged_seurat <- merged_seurat[, !merged_seurat$doublets_consensus.class]

# ==================== 4. 细胞周期评分与校正 ====================

# ==================== 仅使用 DoubletFinder 完成双细胞检测与后续分析 ====================
# 前提：已加载 seurat.data 对象（你的单细胞数据）
# 若未加载，需先加载数据（示例：seurat.data <- readRDS("你的数据路径.rds")）

# ==================== 3.1.3 DoubletFinder 循环处理方案（核心保留） ====================
# 确保样本分组列存在
if (!"sample" %in% colnames(seurat.data@meta.data)) {
  seurat.data$sample <- seurat.data$orig.ident
}

# 获取所有样本名称
samples <- unique(seurat.data$sample)
doublet_results <- list() # 存储每个样本的双细胞预测结果

# 循环处理每个样本
for (samp in samples) {
  cat("\n正在处理样本:", samp, "\n")
  
  # 1. 提取单个样本子集
  sub_seurat <- subset(seurat.data, subset = sample == samp)
  cat("  细胞数:", ncol(sub_seurat), "\n")
  
  # 2. 标准预处理（归一化、找高变基因、标准化、PCA）
  sub_seurat <- NormalizeData(sub_seurat)
  sub_seurat <- FindVariableFeatures(sub_seurat, selection.method = "vst", nfeatures = 2000)
  sub_seurat <- ScaleData(sub_seurat)
  sub_seurat <- RunPCA(sub_seurat)
  
  # 3. 聚类（为 DoubletFinder 提供基础）
  sub_seurat <- FindNeighbors(sub_seurat, dims = 1:20)
  sub_seurat <- FindClusters(sub_seurat, resolution = 0.7)
  
  # 4. 优化 DoubletFinder 参数（pK 自动选择）
  sweep.res <- paramSweep(sub_seurat, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  # 计算预期双细胞数量（基于细胞总数的经验公式）
  nCells_sample <- ncol(sub_seurat)
  nExp_sample <- round(nCells_sample * nCells_sample / 100000 * 0.8) # 每10万细胞对约0.8%双细胞
  
  # 5. 运行 DoubletFinder 检测双细胞
  sub_seurat <- doubletFinder(sub_seurat,
                              PCs = 1:20,
                              pN = 0.25,       # 核心参数：近邻比例（默认0.25即可）
                              pK = optimal_pK, # 自动优化的聚类参数
                              nExp = nExp_sample, # 预期双细胞数量
                              sct = FALSE)     # 若用 SCT 标准化则设为 TRUE
  
  # 6. 提取该样本的预测结果
  df_column_name <- grep("^DF.classifications", colnames(sub_seurat@meta.data), value = TRUE)[1]
  doublet_results[[samp]] <- data.frame(
    Cell = colnames(sub_seurat),
    Doublet_Label = sub_seurat@meta.data[[df_column_name]], # "Doublet" 或 "Singlet"
    Sample = samp,
    row.names = colnames(sub_seurat)
  )
}

# 合并所有样本的结果
library(dplyr)
doublet_df <- bind_rows(doublet_results)
rownames(doublet_df) <- doublet_df$Cell

# 将 DoubletFinder 结果添加到原始 Seurat 对象
seurat.data$Doublet_Label <- "Singlet" # 默认值
matching_cells <- intersect(rownames(doublet_df), colnames(seurat.data))
seurat.data$Doublet_Label[matching_cells] <- doublet_df[matching_cells, "Doublet_Label"]

# ==================== 3.4 双细胞过滤（基于 DoubletFinder 单独结果） ====================
# 无需整合其他方法，直接用 DoubletFinder 结果过滤
# 先将标签转为逻辑值（TRUE=双细胞，FALSE=单细胞），保持与原代码兼容
seurat.data$doublets_consensus.class <- seurat.data$Doublet_Label == "Doublet"

# 输出双细胞检测结果汇总
cat("\n=== DoubletFinder 双细胞预测汇总 ===\n")
print("总体结果:")
print(table(seurat.data$Doublet_Label))
print("\n各样本结果:")
print(table(seurat.data$sample, seurat.data$Doublet_Label))

# 注意：如果之前已对 seurat.data 做过预处理（归一化、高变基因等），可跳过重复步骤
# 若未做过全局预处理，执行以下步骤：
if (!"RNA" %in% Assays(seurat.data) || length(VariableFeatures(seurat.data)) == 0) {
  # 1. 全局归一化（若未做过）
  seurat.data <- NormalizeData(seurat.data)
  # 2. 全局找高变基因（若未做过）
  seurat.data <- FindVariableFeatures(seurat.data, selection.method = "vst", nfeatures = 2000)
  # 3. 全局标准化（若未做过）
  seurat.data <- ScaleData(seurat.data)
}

# 4. 全局 PCA（关键：生成 seurat.data 的 pca 降维结果）
if (!"pca" %in% Reductions(seurat.data)) {
  seurat.data <- RunPCA(seurat.data, features = VariableFeatures(object = seurat.data), seed.use = 1337L)
}

# 可视化双细胞分布（t-SNE）
seurat.data <- RunTSNE(seurat.data, reduction = "pca", dims = 1:20, seed.use = 1337L)
p_doublet <- DimPlot(seurat.data,
                     reduction = "tsne",
                     group.by = "Doublet_Label",
                     cols = c("Singlet" = "lightgrey", "Doublet" = "red"),
                     pt.size = 0.5,
                     order = c("Doublet"), # 双细胞点画在最上层
                     label = FALSE) +
  ggtitle("DoubletFinder 双细胞检测结果") +
  theme(plot.title = element_text(hjust = 0.5))
print(p_doublet)
# 保存图片
ggsave("5.DoubletFinder_Detection.pdf", plot = p_doublet, width = 8, height = 6)
cat("\n双细胞可视化图已保存为：5.DoubletFinder_Detection.pdf\n")

# 过滤双细胞（保留单细胞）
seurat.data_filtered <- seurat.data[, !seurat.data$doublets_consensus.class]
cat("\n过滤后保留细胞数:", ncol(seurat.data_filtered), "\n")
cat("过滤掉的双细胞数:", ncol(seurat.data) - ncol(seurat.data_filtered), "\n")

# ==================== 4. 细胞周期评分与校正（完全保留原流程） ====================
# 4.1 细胞周期基因匹配
# 加载细胞周期基因集（Seurat 内置）
data("cc.genes")

# 分别匹配 S 期和 G2M 期基因（确保与数据中的基因名一致）
g2m.genes <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(seurat.data_filtered))
s.genes <- CaseMatch(search = cc.genes$s.genes, match = rownames(seurat.data_filtered))

# 检查匹配到的基因数量（确保有足够基因用于评分）
cat("\n匹配到的 G2M 期基因数:", length(g2m.genes), "\n")
cat("匹配到的 S 期基因数:", length(s.genes), "\n")

# 4.2 细胞周期评分
seurat.data_filtered <- CellCycleScoring(seurat.data_filtered, 
                                         s.features = s.genes, 
                                         g2m.features = g2m.genes,
                                         assay = "RNA") # 若用 SCT 标准化，改为 assay = "SCT"

# 可视化细胞周期标记基因表达
p_cell_cycle <- RidgePlot(seurat.data_filtered, 
                          features = c("PCNA", "TOP2A", "MCM6", "MKI67"), 
                          ncol = 2) +
  ggtitle("细胞周期标记基因表达分布") +
  theme(plot.title = element_text(hjust = 0.5))
print(p_cell_cycle)
# 保存细胞周期图
ggsave("6.CellCycle_Markers.pdf", plot = p_cell_cycle, width = 10, height = 6)
cat("\n细胞周期可视化图已保存为：6.CellCycle_Markers.pdf\n")

# （可选）细胞周期校正（若后续分析需要去除细胞周期影响）
# seurat.data_filtered <- ScaleData(seurat.data_filtered, vars.to.regress = c("S.Score", "G2M.Score"))

# 保存过滤后的 Seurat 对象（供后续分析使用）
saveRDS(seurat.data_filtered, "seurat_filtered_doublets.rds")
cat("\n过滤后的数据已保存为：seurat_filtered_doublets.rds\n")

# ==================== 4.1 细胞周期基因匹配（仅修改适配部分） ====================
# 1. 加载 Seurat 内置的细胞周期基因集（必须添加，避免 cc.genes 未定义）
if (!exists("cc.genes")) {
  data("cc.genes")
}

# 2. 匹配所有细胞周期基因（S期 + G2M期）- 保留原功能，仅修改对象名
CaseMatch(c(cc.genes$s.genes, cc.genes$g2m.genes), VariableFeatures(seurat.data_filtered))

# 3. 分别匹配 S 期和 G2M 期基因 - 核心修改：merged_seurat → seurat.data_filtered
g2m.genes <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(seurat.data_filtered))
s.genes <- CaseMatch(search = cc.genes$s.genes, match = rownames(seurat.data_filtered))

# （可选但推荐）添加匹配结果检查，方便排查问题
cat("匹配到的 G2M 期基因数:", length(g2m.genes), "\n")
cat("匹配到的 S 期基因数:", length(s.genes), "\n")

# 若匹配基因过少，给出警告（基因名格式可能不兼容）
if (length(g2m.genes) < 5 || length(s.genes) < 5) {
  warning("细胞周期基因匹配数量过少！请检查基因名格式（如大小写、基因符号 vs Ensembl ID）是否一致")
}
# ==================== [F] 4.1 细胞周期基因匹配(原始) ====================

r
# 匹配细胞周期相关基因
CaseMatch(c(cc.genes$s.genes, cc.genes$g2m.genes), VariableFeatures(merged_seurat))

# 分别匹配S期和G2M期基因
g2m.genes <- CaseMatch(search = cc.genes$g2m.genes, match = rownames(merged_seurat))
s.genes <- CaseMatch(search = cc.genes$s.genes, match = rownames(merged_seurat))
# ==================== [F] 4.2 细胞周期评分(原始) ====================

r
# 为每个细胞计算细胞周期评分
merged_seurat <- CellCycleScoring(merged_seurat, s.features = s.genes, g2m.features = g2m.genes)

# 可视化细胞周期相关基因表达
p6_cell_cycle <- RidgePlot(merged_seurat, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
ggsave("545/6.CellCycle_Markers.pdf", plot = p6_cell_cycle, width = 10, height = 6)

图片解释：
•	6.CellCycle_Markers.pdf: 显示细胞周期相关基因(PCNA, TOP2A, MCM6, MKI67)在不同细胞中的表达分布，用于评估细胞周期阶段。
•	若无峰值出现，则忽略细胞周期影响
# ==================== [F] 4.3 细胞周期效应校正（可选）(原始) ====================

r
# 获取所有基因名
all.genes <- rownames(merged_seurat)

# 使用细胞周期评分进行数据校正
merged_seurat <- ScaleData(merged_seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)

# ==================== [F] 5. 环境RNA污染去除（可选）(原始) ====================
# ==================== [F] 5.1 使用decontX去除污染(原始) ====================

r
# 安装并加载celda包
BiocManager::install("celda")
library(celda)

# 获取原始计数数据
counts_data <- merged_seurat@assays$RNA@counts

# 运行decontX去除环境RNA污染
decontX_results <- decontX(counts_data)

# 添加污染分数到元数据
merged_seurat@meta.data$decontX_contamination <- decontX_results$contamination

# 基于污染分数过滤细胞
rt <- merged_seurat@meta.data[merged_seurat@meta.data$decontX_contamination < 0.3, ]

# 更新Seurat对象
merged_seurat <- merged_seurat[, row.names(rt)]

# ==================== [F] 5.2 数据重新标准化(原始) ====================
r
# 重新进行数据标准化和特征选择
merged_seurat <- merged_seurat %>% 
  NormalizeData(verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>%
  ScaleData(verbose = F)

# ==================== 6. 高变基因识别 ====================
# ==================== [F] 6.1 高变基因可视化（原始） ====================

r
# 获取前10个高变基因
top10 <- head(VariableFeatures(seurat.data), 10)

# 绘制高变基因图
plot1 <- VariableFeaturePlot(seurat.data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
p7_variable_features <- plot1 + plot2
ggsave("545/7.Variable_Features.pdf", plot = p7_variable_features, width = 10, height = 6)
print(p7_variable_features)

# ==================== [T] 6.1 高变基因可视化（彻底修复警告 + 窗口问题） ====================
# 1. 绘制高变基因图（用正确参数处理 0 值，避免 log 无限值）
plot1 <- VariableFeaturePlot(
  object = seurat.data,
  assay = "RNA",  # 明确指定 assay（若用 SCT 改为 "SCT"）
  selection.method = "vst",  # 与之前找高变基因的方法一致（确保匹配）
  cols = c("grey", "red"),  # 非高变基因=灰色，高变基因=红色
  pt.size = 0.5  # 点的大小，避免点重叠
)

# 2. 手动替换 x 轴为伪对数转换（解决 0 值导致的无限值警告）
# 原理：将表达量取 log2(表达量 + 1)，既避免 log(0)，又保留表达趋势
plot1 <- plot1 + 
  scale_x_continuous(trans = "pseudo_log", name = "平均表达量 (伪对数转换)") +
  scale_y_continuous(name = "变异系数 (标准化)")

# 3. 获取前10个高变基因（保持原逻辑）
top10 <- head(VariableFeatures(seurat.data), 10)

# 4. 标记高变基因（优化参数，避免标签重叠）
plot2 <- LabelPoints(
  plot = plot1,
  points = top10,
  repel = TRUE,
  xnudge = 0,  # 按提示设置，优化标签位置
  ynudge = 0,
  size = 4,    # 标签文字大小
  fontface = "bold",  # 文字加粗
  color = "black"     # 标签颜色，与红色点对比清晰
)

# 5. 组合图片（确保兼容，避免 patchwork 问题）
library(patchwork)  # 若未安装，先运行：install.packages("patchwork")
p7_variable_features <- plot2  # 直接用带标签的图，无需额外组合

# 6. 保存图片（足够大的尺寸，确保标签完整）
if (!dir.exists("545")) {
  dir.create("545", recursive = TRUE)  # 确保文件夹存在
}
ggsave(
  filename = "545/7.Variable_Features.pdf",
  plot = p7_variable_features,
  width = 12,
  height = 8,
  dpi = 300
)
print(p7_variable_features)
# 7. 强制在独立窗口显示（彻底解决 RStudio 窗口过小问题）
if (interactive()) {
  dev.new(width = 12, height = 8, noRStudioGD = TRUE)  # 新建独立窗口，不依赖 RStudio
}
print(p7_variable_features)

# 8. 输出高变基因名称，方便核对
cat("前10个高变基因：\n")
print(top10)


#图片解释：
#•	7.Variable_Features.pdf: 显示高变基因的识别结果，标记了前10个变异程度最高的基因，这些基因在后续分析中很重要。

# ==================== 7. 降维与批次校正 ====================
# ==================== [T] 7.1 主成分分析(PCA) ====================

r
# 进行PCA降维
seurat.data <- seurat.data %>%
  RunPCA(npcs = 50, verbose = F) %>%
  RunUMAP(reduction = "pca", dims = 1:50, verbose = F)

# 绘制肘部图确定主要成分数量
p8_elbow <- ElbowPlot(seurat.data, ndims = 50)
ggsave("545/8.Elbow_Plot.pdf", plot = p8_elbow, width = 10, height = 10)
print(p8_elbow)
#图片解释：
#•	8.Elbow_Plot.pdf: 肘部图显示每个主成分解释的方差，用于确定在后续分析中使用多少个主成分（选择拐点）。

# ==================== [F] 7.2 批次效应检查（原始） ====================

r
# 检查批次效应
library(patchwork)
options(repr.plot.width = 10, repr.plot.height = 4.5)

p9_batch_effect <- wrap_plots(ncol = 2,
                              DimPlot(merged_seurat, reduction = "pca", group.by = "orig.ident") + NoAxes() + ggtitle("Before_PCA"),
                              DimPlot(merged_seurat, reduction = "umap", group.by = "orig.ident") + NoAxes() + ggtitle("Before_UMAP"),
                              guides = "collect"
)
ggsave("545/9.Batch_Effect_Before.pdf", plot = p9_batch_effect, width = 10, height = 9, dpi = 300)

# ==================== [T] 7.2 批次效应检查（适配 seurat.data 对象） ====================
# 加载 patchwork 包（若未安装：install.packages("patchwork")）
library(patchwork)

# 设置绘图尺寸（可选，不影响核心功能）
options(repr.plot.width = 10, repr.plot.height = 4.5)

# 将 merged_seurat 改为 seurat.data，其他逻辑不变
p9_batch_effect <- wrap_plots(
  ncol = 2,  # 2列布局，同时显示 PCA 和 UMAP
  # PCA 图：按样本（orig.ident）分组，查看批次分布
  DimPlot(
    object = seurat.data, 
    reduction = "pca", 
    group.by = "orig.ident",  # 按原始样本分组，检查是否按批次聚集
    pt.size = 0.3,  # 点大小，避免拥挤
    shuffle = TRUE  # 随机打乱点的绘制顺序，避免重叠遮挡
  ) + NoAxes() + ggtitle("批次效应检查 - PCA"),
  # UMAP 图：同样按样本分组
  DimPlot(
    object = seurat.data, 
    reduction = "umap", 
    group.by = "orig.ident",
    pt.size = 0.3,
    shuffle = TRUE
  ) + NoAxes() + ggtitle("批次效应检查 - UMAP"),
  guides = "collect"  # 收集图例，避免重复
) +
  plot_layout(guides = "collect")  # 统一图例位置，让图更整洁

# 保存图片（确保 545/ 文件夹存在）
if (!dir.exists("545")) {
  dir.create("545", recursive = TRUE)
}
ggsave(
  filename = "545/9.Batch_Effect_Before.pdf",
  plot = p9_batch_effect,
  width = 12,  # 适当增大宽度，容纳图例
  height = 6,
  dpi = 300
)

# 显示图片（用独立窗口，避免尺寸问题）
if (interactive()) {
  dev.new(width = 12, height = 6, noRStudioGD = TRUE)
}
print(p9_batch_effect)

# （可选）补充：统计每个样本的细胞数量，辅助判断批次效应
cat("\n各样本细胞数量统计（批次效应参考）：\n")
print(table(seurat.data$orig.ident))

#图片解释：
#•	9.Batch_Effect_Before.pdf: 显示批次校正前的PCA和UMAP图，不同颜色代表不同样本，用于可视化批次效应。

# ==================== 7.3 Harmony批次校正 ====================

# 设置随机种子保证结果可重复
the.seed <- 1337L
# 使用Harmony进行批次校正
library(harmony)
seurat.data <- seurat.data %>% RunHarmony("orig.ident", plot_convergence = T)

# 基于校正后的数据进行UMAP降维和邻域图构建
n.pcs <- 15
seurat.data <- seurat.data %>%
  RunUMAP(reduction = "harmony", dims = 1:n.pcs, verbose = F, seed.use = the.seed) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n.pcs)

# 可视化批次校正后的效果
p10_after_harmony <- wrap_plots(ncol = 2,
                                DimPlot(seurat.data, reduction = "harmony", group.by = "orig.ident") + NoAxes() + ggtitle("After_PCA (harmony)"),
                                DimPlot(seurat.data, reduction = "umap", group.by = "orig.ident") + NoAxes() + ggtitle("After_UMAP"),
                                guides = "collect"
)
print(p10_after_harmony)
# 比较批次校正前后效果
p11_comparison <- wrap_plots(p9_batch_effect, p10_after_harmony, ncol = 1)
ggsave("545/10.Batch_Correction_Comparison.pdf", plot = p11_comparison, width = 10, height = 9)
ggsave("545/11.After_Harmony.pdf", plot = p10_after_harmony, width = 10, height = 4.5)
print(p11_comparison)
#图片解释：
#•	Batch_Correction_Comparison.pdf: 对比批次校正前后的效果，上方为校正前，下方为校正后。
#•	矫正后细胞混合均匀为矫正有效

# ==================== 8. 细胞聚类分析 ====================
# ==================== 8.1 多分辨率聚类 ====================


# 在不同分辨率下进行聚类
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
  print(res)
  seurat.data <- FindClusters(seurat.data, resolution = res, algorithm = 1) %>% identity()
}

# 绘制聚类树帮助选择合适的分辨率
library("clustree")
p12_clustree <- clustree(seurat.data)
ggsave("545/12.Clustree.pdf", plot = p12_clustree, width = 15, height = 10)
print(p12_clustree)
#   图片解释：
#•	12.Clustree.pdf: 聚类树显示不同分辨率下的聚类关系，帮助选择最合适的分辨率。
#•	一般选择分叉数较少、无杂乱分叉的分辨率

# ==================== 8.2 选择聚类分辨率 ====================

# 选择分辨率0.1进行后续分析
Idents(object = seurat.data) <- "RNA_snn_res.0.1"

# 可视化聚类结果
p13_clusters <- DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.1", label = T) & NoAxes()
ggsave("545/13.Cluster_Resolution_0.1.pdf", plot = p13_clusters, width = 8, height = 7)
print(p13_clusters)
#图片解释：
#•	13.Cluster_Resolution_0.1.pdf: 在UMAP图上显示分辨率0.1下的细胞聚类结果，不同颜色代表不同细胞簇。

# ==================== 9. 细胞类型注释 ====================
# ==================== 9.1 自动注释：SingleR ====================
# 1. 先安装 Bioconductor 核心依赖包（解决 HDF5Array 缺失问题）
# 注意：Bioconductor 包需用 BiocManager 安装，不能用 install.packages
BiocManager::install(c("HDF5Array", "SingleCellExperiment", "SummarizedExperiment"), update = FALSE)

# 2. 验证依赖包是否安装成功
if (requireNamespace("HDF5Array", quietly = TRUE)) {
  cat("HDF5Array 安装成功！\n")
} else {
  stop("HDF5Array 安装失败，请检查网络或 R 版本兼容性！")
}

# 3. 重新安装 celldex（此时依赖已满足）
BiocManager::install("celldex", update = FALSE)

# 4. 验证 celldex 是否安装成功
if (requireNamespace("celldex", quietly = TRUE)) {
  cat("celldex 安装成功！可以正常使用 HumanPrimaryCellAtlasData() 等参考数据集\n")
} else {
  cat("celldex 仍安装失败，建议使用下方的替代方案（跳过 SingleR 自动注释）\n")
}

# 安装并加载SingleR和参考数据库
BiocManager::install("SingleR", update = FALSE)
BiocManager::install("celldex", update = FALSE)
library(SingleR)
library(celldex)
# 1. 安装 SingleR 缺失的依赖包（scrapper + 关联依赖）
# 注意：scrapper 是 Bioconductor 包，必须用 BiocManager 安装
BiocManager::install(c("scrapper", "SingleR", "celldex"), update = FALSE)

# 验证依赖是否安装成功
if (requireNamespace("scrapper", quietly = TRUE)) {
  cat("scrapper 包安装成功！\n")
} else {
  stop("scrapper 安装失败，请检查网络或 R 版本兼容性！")
}

# 加载人类原代细胞图谱参考数据
hpca.se <- celldex::HumanPrimaryCellAtlasData()

# 使用SingleR进行自动细胞类型注释
assay_use <- "RNA"
input_SingleR <- GetAssayData(seurat.data, assay = assay_use, layer = "data")
result_SingleR <- SingleR(test = input_SingleR,
                          ref = hpca.se,
                          labels = hpca.se$label.main,
                          clusters = seurat.data@meta.data$seurat_clusters,
                          quantile = 0.8,
                          fine.tune = F)

# 将注释结果添加到元数据
seurat.data@meta.data$group_cell <- result_SingleR$pruned.labels[match(
  seurat.data@meta.data$seurat_clusters, rownames(result_SingleR))]

# 整理细胞类型信息
celltype <- data.frame(clusterID = rownames(result_SingleR), celltype = result_SingleR$labels)

# 将SingleR注释结果添加到Seurat对象
seurat.data@meta.data$singleR <- "NA"
for(i in 1:nrow(celltype)) {
  seurat.data@meta.data[which(seurat.data$seurat_clusters == celltype$clusterID[i]), 'singleR'] <- celltype$celltype[i]
}

# 可视化SingleR注释结果
p14_singleR <- DimPlot(seurat.data, group.by = "singleR", reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave("545/14.SingleR_Annotation.pdf", plot = p14_singleR, width = 8, height = 7)
print(p14_singleR)


#图片解释：
#•	14.SingleR_Annotation.pdf: 显示基于SingleR自动注释的细胞类型在UMAP图中的分布。

# ==================== 9.2 手动注释：标志基因分析 ====================

# 设置聚类标识
Idents(object = seurat.data) <- "RNA_snn_res.0.1"

# 寻找每个簇的标志基因
allmarkers <- FindAllMarkers(seurat.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 提取每个簇的前5个标志基因
top5_markers <- allmarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# 创建标志基因点图
p15_marker_dotplot <- DotPlot(seurat.data, features = unique(top5_markers$gene)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
  NoLegend() + RotatedAxis()
ggsave("545/15.Marker_Dotplot.pdf", plot = p15_marker_dotplot, width = 30, height = 8)
print(p15_marker_dotplot)
#图片解释：
#•	15.Marker_Dotplot.pdf: 点图显示每个细胞簇中标志基因的表达情况，点的大小表示表达比例，颜色表示平均表达量。

# ==================== 9.3 已知细胞类型标志基因验证 ====================
# 定义已知的细胞类型标志基因
marker1 <- c(
  'AOP4', 'DIO2', 'FPGFR3', 'GFAP', 'ID3', 'LINC00982', 'MT1F', 'SLC14A1', 'SLC1A3', "RYR3", "AQP4",  # 星形胶质细胞
  'COL15A1', 'NOG', 'UNC5B',  # 吊灯样细胞
  'EBF1', 'EMCN', 'ITIH5', 'NOSTRIN', "CLDN5",  # 内皮细胞
  "SLC17A7", "LDB2", "CLSTN2",  # 兴奋性神经元
  "GAD1", "GRIP1", "SYNPR",  # GABA能神经元
  "CDH12", "PAX6",  # 抑制性神经元
  'C3', 'CSF1R', 'CX3CR1', 'TYROBP', "SPP1", "RUNX1", "DOCK8",  # 小胶质细胞
  "GRIN1",  # 神经元
  'KLK6', 'MAG', 'OPALIN', 'PLP1', "ST18", "MBP",  # 少突胶质细胞
  'COL20A1', 'OLIG2', 'PDGFRA', 'PRRX1', "VCAN", "XYLT1"  # 少突胶质前体细胞
)

# 绘制已知标志基因的点图
p16_known_markers <- DotPlot(seurat.data, features = marker1, dot.scale = 8) + RotatedAxis()
ggsave("545/16.Known_Markers.pdf", plot = p16_known_markers, width = 15, height = 12)
print(p16_known_markers)

#图片解释：
#•	16.Known_Markers.pdf: 使用已知的细胞类型标志基因验证聚类结果，帮助确认细胞类型注释的准确性。

# ==================== 9.4 最终细胞类型分配 ====================

# 0.手动检查核心维度
# 维度1：当前Seurat对象的总细胞数（Celltype_Harmony列的行数）
total_cells_meta <- nrow(seurat.data@meta.data)
cat("维度1 - Celltype_Harmony列对应的细胞总数：", total_cells_meta, "\n")

# 维度2：UMAP降维矩阵中的细胞数（有UMAP坐标的细胞数）
total_cells_umap <- nrow(seurat.data@reductions$umap@cell.embeddings)
cat("维度2 - UMAP降维矩阵中的细胞数：", total_cells_umap, "\n")

# 检查差异（报错的核心原因）
cat("两个维度的细胞数差异：", abs(total_cells_meta - total_cells_umap), "\n")

# 额外检查：Celltype_Harmony列的NA值数量（辅助验证）
na_in_annotation <- sum(is.na(seurat.data$Celltype_Harmony))
cat("Celltype_Harmony列的NA值数量：", na_in_annotation, "\n")

# 1. 提取UMAP坐标并筛选出有效坐标的细胞（排除643个NA坐标的细胞）
umap_coords <- Embeddings(seurat.data, reduction = "umap")  # 提取UMAP坐标矩阵
valid_cells <- rownames(umap_coords)[!is.na(umap_coords[,1])]  # 筛选UMAP_1非NA的细胞

# 2. 过滤Seurat对象，只保留有效细胞
seurat.data <- subset(seurat.data, cells = valid_cells)

# 3. 先验证过滤后的细胞数
cat("过滤后实际细胞数：", ncol(seurat.data), "\n")

# 4. 注释列初始化：直接设为Unknown，彻底杜绝NA（核心修改）
seurat.data$Celltype_Harmony <- "Unknown"

# 5. 提取过滤后的聚类ID（避免用旧的Idents匹配）
idents <- as.character(Idents(seurat.data))

# 6. 给已知聚类赋值（未提及的聚类默认Unknown，无NA）
seurat.data$Celltype_Harmony[idents %in% c("0")] <- "Oligodendrocyte"  # 少突胶质细胞
seurat.data$Celltype_Harmony[idents %in% c("1")] <- "Astrocytes"       # 星形胶质细胞
seurat.data$Celltype_Harmony[idents %in% c("4")] <- "OPC"              # 少突胶质前体细胞
seurat.data$Celltype_Harmony[idents %in% c("5")] <- "Microglial cell"  # 小胶质细胞
seurat.data$Celltype_Harmony[idents %in% c("2", "3")] <- "Excitatory neuron"  # 兴奋性神经元
seurat.data$Celltype_Harmony[idents %in% c("6", "7", "8")] <- "Inhibitory neuron"  # 抑制性神经元
# 9/10/11/12/13/14等聚类已默认是Unknown，无需单独赋值

# 7. 验证：注释列无NA值（必须为0）
cat("注释列NA值数量：", sum(is.na(seurat.data$Celltype_Harmony)), "\n") # 必须显示0

# 可视化最终细胞类型注释结果
p17_final_annotation <- DimPlot(seurat.data, group.by = "Celltype_Harmony", reduction = "umap", label = TRUE)
ggsave("545/17.Final_Celltype_Annotation.pdf", plot = p17_final_annotation, width = 8, height = 7)
print(p17_final_annotation)
#图片解释：
#•	17.Final_Celltype_Annotation.pdf: 显示最终手动注释的细胞类型在UMAP图中的分布，包含主要的脑细胞类型。

# ==================== 10. 差异表达分析 ====================

# 10.1 疾病状态差异分析

# 进行AD组与HC组之间的差异基因分析
# FindMarkers默认基于 Seurat 对象的Idents（细胞身份）分组，而非直接读取meta.data中的condition列
# 所以哪怕seurat.data预处理部分映射无误，这里也必然会报错，我不知道这段为什么这么写
all_cell_deg <- FindMarkers(
  object = seurat.data,        # Seurat对象
  ident.1 = "AD",                # 实验组：阿尔茨海默病
  ident.2 = "HC",                # 对照组：健康
  only.pos = FALSE,              # 保留上调和下调基因
  logfc.threshold = 0.1,         # log2FC阈值
  min.pct = 0.25,                # 基因在组内的最低表达比例
  assay = "RNA"                  # 使用的assay
)

# 处理方法：
# 1. 查看当前Idents的类型（是聚类ID还是condition）
cat("当前Idents的取值：", paste(unique(Idents(seurat.data)), collapse = ", "), "\n")
# 若输出是0、1、2...（聚类ID），就是核心问题；若输出AD/HC，则跳过第二步

# 2. 验证condition列确实存在且有值
cat("condition列是否存在：", "condition" %in% colnames(seurat.data@meta.data), "\n")
cat("condition列的取值：", paste(unique(seurat.data$condition), collapse = ", "), "\n")
# 正常输出应为：TRUE + AD, HC

# 切换Idents到condition列
Idents(seurat.data) <- "condition"

# 验证切换结果
cat("切换后Idents的取值：", paste(unique(Idents(seurat.data)), collapse = ", "), "\n")
# 必须输出：AD, HC

# 重新运行差异基因分析（此时可正常识别AD/HC）
all_cell_deg <- FindMarkers(
  object = seurat.data,        
  ident.1 = "AD",                # 实验组：阿尔茨海默病
  ident.2 = "HC",                # 对照组：健康
  only.pos = FALSE,              # 保留上调和下调基因
  logfc.threshold = 0.1,         # log2FC阈值
  min.pct = 0.25,                # 基因在组内的最低表达比例
  assay = "RNA"                  # 使用的assay
)

# 验证结果
head(all_cell_deg)
cat("差异基因总数：", nrow(all_cell_deg), "\n")


# 设置差异表达阈值
logFC_threshold <- 0.10
pvalue_threshold <- 0.05

# 标记差异表达基因
deg <- all_cell_deg
deg$change <- ifelse(deg$p_val_adj < pvalue_threshold & deg$avg_log2FC > logFC_threshold, "up",
                     ifelse(deg$p_val_adj < pvalue_threshold & deg$avg_log2FC < -logFC_threshold, "down", "stable"))

# 统计差异基因数量
table(deg$change)
# down stable     up 
# 999    490    505 
# up:上调基因：这些基因在 AD 组细胞中的表达量，显著高于 HC 组（log2FC>0.1，校正 P<0.05）
# down：下调基因：这些基因在 AD 组细胞中的表达量，显著低于 HC 组（log2FC<-0.1，校正 P<0.05）
# stable：无差异基因：这些基因在 AD 和 HC 组中的表达量无统计学显著差异（要么倍数不够，要么 P 值不显著）


