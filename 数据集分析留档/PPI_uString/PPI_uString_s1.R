# ========================= 第一步：安装并加载所需R包 =========================
# 安装核心包（首次运行需安装，已安装则跳过）
if (!require("STRINGdb", quietly = TRUE)) {
  install.packages("STRINGdb")  # String数据库R接口
}
if (!require("igraph", quietly = TRUE)) {
  install.packages("igraph")    # 网络分析核心包
}
if (!require("ggraph", quietly = TRUE)) {
  install.packages("ggraph")    # 高颜值网络可视化
}
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")   # 可视化基础包
}
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")     # 数据处理
}
if (!require("tidyr", quietly = TRUE)) {
  install.packages("tidyr")     # 数据格式转换
}

# 加载包
library(STRINGdb)
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(tidyr)

# ========================= 第二步：准备上调/下调基因列表 =========================
# PPI分析核心输入：上调基因列表 + 下调基因列表（需为HGNC SYMBOL格式）
# ---------------------- 2.1 模拟示例基因列表（用户替换为自己的真实数据） ----------------------
set.seed(123)  # 固定随机种子，结果可重复
# 模拟上调基因（20个，示例）
up_genes <- c("TP53", "MYC", "EGFR", "KRAS", "AKT1", "CDKN1A", "MDM2", "VEGFA", 
              "CCND1", "MAPK1", "BCL2", "PTEN", "PIK3CA", "RB1", "TNF", "IL6", 
              "CASP3", "FOS", "JUN", "STAT3")
# 模拟下调基因（15个，示例）
down_genes <- c("CD44", "MET", "NOTCH1", "SMAD4", "BRCA1", "ATM", "CHEK2", "PTGS2", 
                "HIF1A", "MMP9", "CXCL8", "ERK2", "AKT2", "MYB", "FOXO1")

# ---------------------- 2.2 真实基因列表加载（用户替换核心代码） ----------------------
# 方式1：从CSV文件加载（示例：第一列是基因名，列名为"Gene"）
# up_genes <- read.csv("up_genes.csv")$Gene    # 上调基因列表
# down_genes <- read.csv("down_genes.csv")$Gene  # 下调基因列表
# 方式2：从TXT文件加载（每行一个基因名）
# up_genes <- read.table("up_genes.txt")$V1
# down_genes <- read.table("down_genes.txt")$V1

# ---------------------- 2.3 基因列表预处理（关键，避免报错） ----------------------
# 合并上调/下调基因（用于构建整体PPI网络）
all_genes <- c(up_genes, down_genes)
# 去重（避免重复基因干扰分析）
all_genes <- unique(all_genes)
# 转换为字符型（确保格式正确）
all_genes <- as.character(all_genes)
# 过滤空值/NA（避免连接数据库报错）
all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]

# ========================= 第三步：连接String数据库并获取PPI网络数据 =========================
# ---------------------- 3.1 初始化String数据库连接 ----------------------
# string_db <- STRINGdb$new(
#   version = "11.5",          # String数据库版本（最新为11.5）
#   species = 9606,            # 物种编号：人类=9606，小鼠=10090，大鼠=10116
#   score_threshold = 400,     # 相互作用置信度阈值（0-1000，越高越严格，400=中等置信度）
#   input_directory = getwd()  # 数据缓存目录（当前工作目录）
# )

# 初始化String数据库连接（简化版，自动加载最新版本）
string_db <- STRINGdb$new(
  species = 9606,            # 物种：人类（9606），小鼠替换为10090
  score_threshold = 400,     # 置信度阈值：400（中等），建议范围400-700
  input_directory = getwd()  # 缓存目录（无需修改）
)

# ---------------------- 3.2 将基因列表映射到String数据库 ----------------------
# 映射基因名（解决基因名别名问题，确保匹配String数据库）
mapped_genes <- string_db$map(data.frame(Gene = all_genes), 
                              "Gene",  # 基因名列名
                              removeUnmappedRows = TRUE)  # 移除无法映射的基因
# 提取成功映射的基因（用于后续分析）
mapped_gene_list <- unique(mapped_genes$STRING_id)
cat("成功映射到String数据库的基因数：", length(mapped_gene_list), "/", length(all_genes), "\n")

# ---------------------- 3.3 获取PPI相互作用网络数据 ----------------------
# 调取基因间的相互作用关系（edges）
ppi_edges <- string_db$get_interactions(mapped_gene_list)
# 查看PPI数据结构：from, to, combined_score（相互作用置信度）
print(head(ppi_edges))

# ---------------------- 3.4 过滤低置信度相互作用（可选但推荐） ----------------------
# 只保留combined_score ≥ 400的相互作用（与初始化阈值一致）
ppi_edges_filtered <- ppi_edges[ppi_edges$combined_score >= 400, ]
# 移除重复的相互作用（A-B和B-A视为同一对，只保留一条）
ppi_edges_filtered <- ppi_edges_filtered %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(from, to)), collapse = "_")) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

# ========================= 第四步：构建PPI网络并标记上调/下调基因 =========================
# ---------------------- 4.1 构建igraph网络对象 ----------------------
# 从edges构建无向网络
ppi_network <- graph_from_data_frame(
  d = ppi_edges_filtered[, c("from", "to")],  # 边列表（from=源节点，to=目标节点）
  directed = FALSE,                            # 无向网络（蛋白相互作用无方向）
  vertices = data.frame(id = mapped_gene_list)  # 节点列表（映射后的基因）
)

# ---------------------- 4.2 为节点添加属性：上调/下调/无差异 ----------------------
# 先将String ID转换回基因名（便于标注）
node_gene_map <- mapped_genes %>%
  select(STRING_id, Gene) %>%
  distinct(STRING_id, .keep_all = TRUE)
# 为网络节点添加基因名属性
V(ppi_network)$gene_name <- node_gene_map$Gene[match(V(ppi_network)$name, node_gene_map$STRING_id)]
# 为节点添加"表达趋势"属性（up/down/none）
V(ppi_network)$expression <- ifelse(V(ppi_network)$gene_name %in% up_genes, 
                                    "Up", 
                                    ifelse(V(ppi_network)$gene_name %in% down_genes, 
                                           "Down", "None"))
# 过滤掉"None"节点（仅保留上调/下调基因）
ppi_network_filtered <- induced_subgraph(ppi_network, 
                                         V(ppi_network)$expression %in% c("Up", "Down"))

# ---------------------- 4.3 计算节点核心性指标（筛选核心基因） ----------------------
# 计算度中心性（Degree）：节点连接的边数，越高说明该基因在网络中越核心
V(ppi_network_filtered)$degree <- degree(ppi_network_filtered)
# 计算介数中心性（Betweenness）：节点在网络中作为"桥梁"的程度，越高越关键
V(ppi_network_filtered)$betweenness <- betweenness(ppi_network_filtered)
# 计算紧密中心性（Closeness）：节点到其他节点的平均最短路径，越高越核心
V(ppi_network_filtered)$closeness <- closeness(ppi_network_filtered)

# 提取核心基因（度中心性前10%的基因）
core_genes <- V(ppi_network_filtered)$gene_name[V(ppi_network_filtered)$degree >= quantile(V(ppi_network_filtered)$degree, 0.9)]
cat("PPI网络核心基因（度中心性前10%）：", paste(core_genes, collapse = ", "), "\n")

# ========================= 第五步：PPI网络可视化（核心） =========================
# ---------------------- 5.1 基础PPI网络可视化（区分上调/下调） ----------------------
# 设置可视化主题
set.seed(123)  # 固定布局，确保图的位置不变
ggraph(ppi_network_filtered, layout = "fr") +  # "fr"布局：力导向布局（美观）
  # 绘制边（相互作用）：灰色，半透明
  geom_edge_link(color = "gray50", alpha = 0.5, width = 0.5) +
  # 绘制节点（基因）：按表达趋势着色，大小按度中心性（核心基因更大）
  geom_node_point(aes(color = expression, size = degree), alpha = 0.8) +
  # 标记核心基因名称（度中心性前10%）
  geom_node_text(aes(label = ifelse(degree >= quantile(degree, 0.9), gene_name, "")), 
                 size = 3, repel = TRUE) +  # repel避免文字重叠
  # 设置颜色：上调=红色，下调=蓝色
  scale_color_manual(values = c("Up" = "#e74c3c", "Down" = "#3498db")) +
  # 设置节点大小范围
  scale_size(range = c(2, 8)) +
  # 标题和标签
  labs(title = "PPI Network of Up/Down Regulated Genes",
       color = "Expression Trend",
       size = "Degree Centrality") +
  # 美化主题
  theme_graph() +  # 网络专用主题（无背景/坐标轴）
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right")

# 保存网络图（高清PDF）
ggsave("PPI_Network_Up_Down_Genes.pdf", 
       width = 12, height = 10, dpi = 300)

# ---------------------- 5.2 核心基因子网络可视化（可选） ----------------------
# 提取核心基因的子网络（仅展示度前10%的基因及相互作用）
core_gene_nodes <- V(ppi_network_filtered)$name[V(ppi_network_filtered)$gene_name %in% core_genes]
core_subnetwork <- induced_subgraph(ppi_network_filtered, core_gene_nodes)

# 绘制核心子网络
set.seed(123)
ggraph(core_subnetwork, layout = "fr") +
  geom_edge_link(color = "gray50", alpha = 0.7, width = 1) +
  geom_node_point(aes(color = expression, size = degree), alpha = 0.9) +
  geom_node_text(aes(label = gene_name), size = 4, repel = TRUE) +
  scale_color_manual(values = c("Up" = "#e74c3c", "Down" = "#3498db")) +
  scale_size(range = c(5, 10)) +
  labs(title = "Core Genes Subnetwork (Top 10% Degree)",
       color = "Expression Trend",
       size = "Degree Centrality") +
  theme_graph() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# 保存核心子网络
ggsave("Core_Genes_PPI_Subnetwork.pdf", 
       width = 10, height = 8, dpi = 300)

# ========================= 第六步：保存分析结果 =========================
# ---------------------- 6.1 保存PPI网络数据 ----------------------
# 提取节点属性（基因名、表达趋势、核心性指标）
node_attributes <- data.frame(
  Gene_Name = V(ppi_network_filtered)$gene_name,
  STRING_ID = V(ppi_network_filtered)$name,
  Expression_Trend = V(ppi_network_filtered)$expression,
  Degree_Centrality = V(ppi_network_filtered)$degree,
  Betweenness_Centrality = V(ppi_network_filtered)$betweenness,
  Closeness_Centrality = V(ppi_network_filtered)$closeness
)
# 保存节点属性
write.csv(node_attributes, "PPI_Node_Attributes.csv", row.names = FALSE)

# ---------------------- 6.2 保存边（相互作用）数据 ----------------------
edge_attributes <- as_data_frame(ppi_network_filtered, what = "edges") %>%
  left_join(node_attributes[, c("STRING_ID", "Gene_Name")], by = c("from" = "STRING_ID")) %>%
  rename(Gene_From = Gene_Name) %>%
  left_join(node_attributes[, c("STRING_ID", "Gene_Name")], by = c("to" = "STRING_ID")) %>%
  rename(Gene_To = Gene_Name) %>%
  select(Gene_From, Gene_To, everything())
# 保存边数据
write.csv(edge_attributes, "PPI_Edge_Attributes.csv", row.names = FALSE)

# ---------------------- 6.3 保存核心基因列表 ----------------------
core_gene_df <- data.frame(Core_Gene = core_genes,
                           Degree = V(ppi_network_filtered)$degree[match(core_genes, V(ppi_network_filtered)$gene_name)])
write.csv(core_gene_df, "Core_Genes_List.csv", row.names = FALSE)

#专注于蛋白相互作用的公共数据库，STRINGdb包是官方 R 接口，支持直接调取相互作用数据。
#物种切换：
#人类：species = 9606
#小鼠：species = 10090
#大鼠：species = 10116
#需根据基因列表物种对应修改，否则基因映射失败。
#置信度阈值：
#score_threshold = 400（中等置信度，适合初步分析）
#严格分析可设700（高置信度，仅保留强相互作用）
#宽松分析可设200（低置信度，包含更多弱相互作用）
#结果：
#节点颜色：红色 = 上调基因，蓝色 = 下调基因；
#节点大小：越大表示度中心性越高（网络核心基因）；
#核心基因：度 / 介数 / 紧密中心性高的基因，通常是调控通路的关键节点。

# 从 R 中导出 Cytoscape 专属的 PPI 网络文件
# 加载之前的PPI结果
edge_df <- read.csv("PPI_Edge_Attributes.csv")  # 边列表
node_df <- read.csv("PPI_Node_Attributes.csv")  # 节点列表

# 提取边列表核心列（重命名为Cytoscape默认识别的列名）
cytoscape_edge <- edge_df[, c("Gene_From", "Gene_To")]
# 提取节点列表核心列（保留基因名+表达趋势，用于着色）
cytoscape_node <- node_df[, c("Gene_Name", "Expression_Trend")]

# 保存为Cytoscape专属文件（UTF-8格式，避免中文乱码，无行名）
write.csv(cytoscape_edge, "Cytoscape_PPI_Edge.csv", row.names = F, fileEncoding = "UTF-8")
write.csv(cytoscape_node, "Cytoscape_PPI_Node.csv", row.names = F, fileEncoding = "UTF-8")