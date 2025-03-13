
############################################### 富集分析函数 #############################################
# GO 富集分析
GOEnrichment <- function(geneDf, pAdjustMethod = "BH", pvalueCutoff = 1, 
                         qvalueCutoff = 1, topNumber = 20, csvBaseName = NULL) {
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(tidyverse)
  
  ego <- enrichGO(gene = geneDf$SYMBOL, 
                  OrgDb = org.Mm.eg.db, 
                  keyType = "SYMBOL",
                  ont = "all", 
                  pAdjustMethod = pAdjustMethod,
                  pvalueCutoff = pvalueCutoff, 
                  qvalueCutoff = qvalueCutoff)
  
  go <- as.data.frame(ego)
  go$Term <- paste(go$ID, go$Description, sep = ': ')
  go <- go[order(go$ONTOLOGY, go$p.adjust, method = "radix", decreasing = c(TRUE, FALSE)), ]
  go$Term <- factor(go$Term, levels = go$Term)
  
  # Top term of each group
  top_go <- go %>%
    group_by(ONTOLOGY) %>%
    slice_head(n = topNumber) %>%
    ungroup()
  
  write.csv(go, paste0("GOAll", csvBaseName, ".csv"), row.names = FALSE)
  write.csv(top_go, paste0("GOTop", topNumber, csvBaseName, ".csv"), row.names = FALSE)
  
  return(list(go = go, top_go = top_go))
}

KEGGEnrichment <- function(geneDf, pAdjustMethod = "BH", pvalueCutoff = 1, 
                           qvalueCutoff = 1, topNumber = 30, csvBaseName = NULL) {
  
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(tidyverse)
  
  # 检查 geneDf 是否包含 ENTREZID 列
  if (!"ENTREZID" %in% colnames(geneDf)) {
    if ("SYMBOL" %in% colnames(geneDf)) {
      convert <- bitr(geneDf$SYMBOL,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)
      
      # 合并转换后的 ENTREZID 到 geneDf
      geneDf <- merge(geneDf, convert, by.x = "SYMBOL", by.y = "SYMBOL", all.x = TRUE)
    } else {
      stop("输入数据框中需要含有 ENTREZID 列或 SYMBOL 列。")
    }
  }
  
  # 进行 KEGG 富集分析
  ekegg <- enrichKEGG(gene = geneDf$ENTREZID, 
                      organism = "mmu", 
                      pAdjustMethod = pAdjustMethod,
                      pvalueCutoff = pvalueCutoff, 
                      qvalueCutoff = qvalueCutoff)
  
  # 处理结果
  ekegg@result$Pathway <- paste(ekegg@result$ID, ekegg@result$Description, sep = ': ')
  ekegg2 <- setReadable(ekegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  
  # 提取前 topNumber 个通路，确保不会生成 NA 行
  if (nrow(ekegg2@result) >= topNumber) {
    top_kegg <- ekegg2@result[1:topNumber, ]
  } else {
    top_kegg <- ekegg2@result  # 如果行数不足，提取所有行
  }
  
  # 保存结果到 CSV 文件
  write.csv(ekegg2, paste0("KEGGAll", csvBaseName, ".csv"), row.names = F)
  write.csv(top_kegg, paste0("KEGGTop", topNumber, csvBaseName, ".csv"), row.names = F)
  
  # 返回结果
  return(list(ekegg = ekegg2, top_kegg = top_kegg))
}



# 绘图
plotGO <- function(ego, color = "Common", seed = NULL, facetGrid = TRUE, discreteWidth = 150, 
                   title = "GO Enrichment Analysis", pFileBaseName = NULL) {
  library(clusterProfiler)
  library(ggplot2)
  library(tidyverse)
  
  if (color == "Common") {
    colors <- c("#8DA1CB", "#FD8D62", "#66C3A5")
  } else if (color == "Random") {
    colors <- randomColor(3, seed)
    print(paste("Use colors:", paste(colors, collapse = " ")))
  }
  
  p <- ggplot(ego, aes(Term, Count)) + 
    geom_col(aes(fill = ONTOLOGY), width = 0.8, show.legend = TRUE) + 
    scale_fill_manual(values = colors) + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = discreteWidth)) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
    coord_flip() + 
    theme_bw() + 
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
          axis.text = element_text(color = "black", size = 12), legend.text = element_text(size = 13)) + 
    labs(x = 'GO Terms', title = title)
  
  # 分面
  if (facetGrid == "TRUE") {
    p <- p + facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')
  }
  
  return(p)
}

plotKEGG <- function(ekegg, color = "Common", seed = NULL, facetGrid = TRUE, discreteWidth = 100, 
                     title = "KEGG Enrichment Analysis", pFileBaseName = NULL) {
  library(clusterProfiler)
  library(ggplot2)
  library(tidyverse)
  
  if (color == "Common") {
    colors <- c("#8DA1CB", "#FD8D62", "#66C3A5")
  } else if (color == "Random") {
    colors <- randomColor(2, seed)
    print(paste("Use colors:", paste(colors, collapse = " ")))
  }
  
  # p <- barplot(ekegg) + 
  #   labs(y = "KEGG Pathways", title = title) + 
  #   scale_y_discrete(labels = function(x) str_wrap(x, width = discreteWidth)) + 
  #   scale_fill_manual(colors)
  #   # scale_fill_gradient(low = colors[1], high = colors[2])
  
  ekegg$Pathway <- factor(ekegg$Pathway, levels = ekegg$Pathway[order(ekegg$p.adjust)])
  
  p <- ggplot(ekegg, aes(Count, Pathway)) +
    geom_col(aes(fill = p.adjust), width = 0.8, show.legend = TRUE) +
    scale_fill_gradient(low = colors[1], high = colors[2]) +
    labs(y = "KEGG Pathways", title = title) +
    theme_bw() +
    theme(axis.text = element_text(color = "black", size = 12))
  
  return(p)
}
