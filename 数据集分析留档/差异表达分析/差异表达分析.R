setwd("E:/Drug_Design/DEGs")
library(tidyverse)
library(limma)


exp <- read.csv("series_matrix.csv", row.names = 1)

exp <- exp[rowMeans(exp) > 1, ]

if(max(exp, na.rm = TRUE) > 15) {
  exp <- log2(exp + 1)
  cat("log2finish！！！\n")
}

boxplot(exp, outline = FALSE, notch = FALSE, las = 2, main = "before")
exp <- normalizeBetweenArrays(exp)
boxplot(exp, outline = FALSE, notch = FALSE, las = 2, main = "after")

# 基因名转换
annotation <- read.csv("GPL.csv", stringsAsFactors = FALSE)

id_to_symbol <- setNames(annotation$Gene.Symbol, annotation$ID)

exp <- exp[rownames(exp) %in% annotation$ID, ] %>%
  as.matrix() %>%
  {
    # 获取匹配的基因符号
    symbols <- id_to_symbol[rownames(.)]
    
    # 处理多基因注释（保留第一个基因名）
    symbols_processed <- sapply(strsplit(symbols, " /// "), function(x) {
      if(length(x) > 0 && !is.na(x[1])) x[1] else NA
    })
    
    # 重复基因符号取平均值
    aggregated <- aggregate(., by = list(GeneSymbol = symbols_processed), FUN = mean, na.rm = TRUE)
    rownames(aggregated) <- aggregated$GeneSymbol
    aggregated$GeneSymbol <- NULL
    
    # 移除无效行
    valid_rows <- !is.na(rownames(aggregated)) & rownames(aggregated) != "" & !grepl("^\\s*$", rownames(aggregated))
    as.data.frame(aggregated[valid_rows, ])
  }

write.csv(exp, "new_series_matrix.csv")

#差异分析------------------------------

fen<-read.csv("group.csv",row.names = 1)#读入表达矩阵表头的分类文件

group_list <- factor(fen$group,levels = c("Control","Test"))
design <- model.matrix(~group_list)
con <- lmFit(exp,design)
con2 <- eBayes(con)
DEG1 <- topTable(con2, coef = 2, number = Inf)
DEG2 = na.omit(DEG1) 

logFC_cut = 0.25
p_cut = 0.05

type1 = (DEG2$adj.P.Val < p_cut)&(DEG2$logFC < -logFC_cut)
type2 = (DEG2$adj.P.Val < p_cut)&(DEG2$logFC > logFC_cut)

DEG2$type = ifelse(type1,"Down",ifelse(type2,"Up","NOT"))

head(DEG2)

write.csv(DEG2,"DEGs.csv")

# 火山图绘制---------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(tibble) 

volcano_data <- DEG2 %>%
  rownames_to_column("gene") %>%  
  mutate(
    Regulation = case_when(
      logFC > 0.25 & adj.P.Val < 0.05 ~ "Up",
      logFC < -0.25 & adj.P.Val < 0.05 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

sig_up <- sum(volcano_data$Regulation == "Up")
sig_down <- sum(volcano_data$Regulation == "Down")

volcano_data$Regulation <- factor(volcano_data$Regulation, 
                                  levels = c("Up", "Down", "Not significant"))

color_values <- c("Up" = "#E9687A", "Down" = "#B6B3D6", "Not significant" = "#D5D1D1")
legend_labels <- c(
  sprintf("Up (n = %d)", sig_up),
  sprintf("Down (n = %d)", sig_down),
  "Not significant"
)

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = Regulation)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(
    values = color_values,
    labels = legend_labels
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  
  geom_text_repel(
    data = volcano_data %>% 
      filter(Regulation != "Not significant") %>%
      arrange(adj.P.Val) %>%
      group_by(Regulation) %>%
      slice_head(n = 10),
    aes(label = gene),  
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    segment.color = "grey50"
  ) +
  
  labs(
    title = "Differential Gene Expression Analysis",
    subtitle = "|log2FC| > 0.25, Adjusted p-value < 0.05",
    x = "log2 Fold Change (R vs NR)", 
    y = "-log10(Adjusted p-value)",
    color = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_x_continuous(breaks = seq(floor(min(volcano_data$logFC)), 
                                  ceiling(max(volcano_data$logFC)), 
                                  by = 0.5))

ggsave("Gene_Expression_Volcano.pdf", width = 10, height = 8, dpi = 300)

