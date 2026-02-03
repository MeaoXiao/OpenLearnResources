# 安装（首次使用运行）
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"))

# 加载
library(clusterProfiler)
library(org.Hs.eg.db)  # 人类注释库（小鼠换org.Mm.eg.db）
library(enrichplot)
library(ggplot2)
diff_result <- DEG2

# 1. 定义筛选阈值（可根据需求调整）
fc_cutoff <- 1    # log2差异倍数阈值（|log2FC|>1）
p_cutoff <- 0.1  # 校正后P值阈值

# 2. 筛选上调基因（logFC>1，校正P值<0.05）
up_genes <- rownames(diff_result[diff_result$logFC > fc_cutoff & diff_result$adj.P.Val < p_cutoff, ])
# 3. 筛选下调基因（logFC<-1，校正P值<0.05）
down_genes <- rownames(diff_result[diff_result$logFC < -fc_cutoff & diff_result$adj.P.Val < p_cutoff, ])

# 检查数量（至少5个基因才有富集意义）
cat("上调基因数量：", length(up_genes), "\n")
cat("下调基因数量：", length(down_genes), "\n")

# ========== 上调基因GO分析 ==========
go_up <- enrichGO(
  gene = up_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",  # 直接用Gene Symbol，无需ID转换
  ont = "ALL",         # 分析BP/CC/MF全部3类（也可指定BP）
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.05,
  readable = TRUE      # 结果中显示基因名，更易读
)
go_up_result <- as.data.frame(go_up)  # 转为数据框

# ========== 下调基因GO分析 ==========
go_down <- enrichGO(
  gene = down_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.05,
  readable = TRUE
)
go_down_result <- as.data.frame(go_down)

# 查看结果（可选）
cat("上调基因GO富集条目数：", nrow(go_up_result), "\n")
cat("下调基因GO富集条目数：", nrow(go_down_result), "\n")

# 先将Gene Symbol转为Entrez ID（KEGG必需）
# 上调基因ID转换
up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
# 下调基因ID转换
down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

# ========== 上调基因KEGG分析 ==========
kegg_up <- enrichKEGG(
  gene = up_entrez,
  organism = "hsa",  # 人类（小鼠换mmu）
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr"
)
kegg_up_result <- as.data.frame(kegg_up)

# ========== 下调基因KEGG分析 ==========
kegg_down <- enrichKEGG(
  gene = down_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr"
)
kegg_down_result <- as.data.frame(kegg_down)

# 查看结果（可选）
cat("上调基因KEGG富集条目数：", nrow(kegg_up_result), "\n")
cat("下调基因KEGG富集条目数：", nrow(kegg_down_result), "\n")

# 1. GO富集气泡图（上调+下调分开画）
# 上调GO
dotplot(go_up, showCategory = 10, title = "Up-regulated Genes - GO Enrichment")
ggsave("GO_up_bubble.png", width = 10, height = 8)

# 下调GO
dotplot(go_down, showCategory = 10, title = "Down-regulated Genes - GO Enrichment")
ggsave("GO_down_bubble.png", width = 10, height = 8)

# 2. KEGG富集气泡图（上调+下调分开画）
# 上调KEGG
dotplot(kegg_up, showCategory = 10, title = "Up-regulated Genes - KEGG Enrichment")
ggsave("KEGG_up_bubble.png", width = 10, height = 8)

# 下调KEGG
dotplot(kegg_down, showCategory = 10, title = "Down-regulated Genes - KEGG Enrichment")
ggsave("KEGG_down_bubble.png", width = 10, height = 8)

# 保存GO结果
write.csv(go_up_result, "GO_up_regulated.csv", row.names = FALSE)
write.csv(go_down_result, "GO_down_regulated.csv", row.names = FALSE)

# 保存KEGG结果
write.csv(kegg_up_result, "KEGG_up_regulated.csv", row.names = FALSE)
write.csv(kegg_down_result, "KEGG_down_regulated.csv", row.names = FALSE)