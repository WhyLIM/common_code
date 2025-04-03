
############################################### 富集分析函数 #############################################
#' GO富集分析函数
#'
#' 对基因列表进行GO富集分析，支持人类和小鼠数据，结果自动保存为 CSV 文件。
#'
#' @param geneDf 数据框，必须包含 SYMBOL 列
#' @param pAdjustMethod p 值校正方法，默认为"BH"（Benjamini-Hochberg）
#' @param pvalueCutoff p 值阈值，默认为 1（不过滤）
#' @param qvalueCutoff q 值阈值，默认为 1（不过滤）
#' @param topNumber 每组（BP/MF/CC）提取条目数，默认为 20
#' @param suffix 输出文件名的后缀标识（自动添加下划线连接），默认为 NULL
#' @param species 物种选择，可选 "human" 或 "mouse"，默认为"human"
#' @param outputDir 输出文件保存目录，默认为当前工作目录（"."）
#'
#' @return 返回包含两个元素的列表：
#' \itemize{
#'   \item go - 完整 GO 富集结果数据框
#'   \item top_go - 每组（BP/MF/CC）的 top 条目结果
#' }
#'
#' @examples
#' \dontrun{
#' # 人类基因分析（保存到当前目录）
#' GOEnrichment(human_geneDf, suffix = "patient_samples")
#'
#' # 小鼠基因分析（保存到指定目录）
#' GOEnrichment(mouse_geneDf, 
#'              species = "mouse", 
#'              outputDir = "./GO")
#'
#' # 保存返回变量
#' result <- GOEnrichment(geneDf, suffix = NULL)
#' }
#'
#' @export
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom tidyverse %>%
GOEnrichment <- function(geneDf, 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 1, 
                         qvalueCutoff = 1, 
                         topNumber = 20, 
                         suffix = NULL, 
                         species = c("human", "mouse"),
                         outputDir = ".") {
  
  # 隐藏包加载时自动显示的启动消息
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(tidyverse)
    # 根据物种加载对应的OrgDb
    species <- match.arg(species)
    if (species == "mouse") {
      library(org.Mm.eg.db)
      OrgDb <- org.Mm.eg.db
    } else {
      library(org.Hs.eg.db)
      OrgDb <- org.Hs.eg.db
    }
  })
  
  # 检查输出目录是否存在，不存在则创建
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
    message("创建输出目录: ", outputDir)
  }
  
  # 参数检查
  if (!"SYMBOL" %in% colnames(geneDf)) {
    stop("输入数据框必须包含 SYMBOL 列")
  }
  
  ego <- enrichGO(gene = geneDf$SYMBOL, 
                  OrgDb = OrgDb, 
                  keyType = "SYMBOL",
                  ont = "all", 
                  pAdjustMethod = pAdjustMethod,
                  pvalueCutoff = pvalueCutoff, 
                  qvalueCutoff = qvalueCutoff)
  
  go <- as.data.frame(ego)
  go$Term <- paste(go$ID, go$Description, sep = ': ')
  # 按 ONTOLOGY 分组后按 p.adjust 升序排序
  go <- go[order(go$ONTOLOGY, go$p.adjust, method = "radix", decreasing = c(TRUE, FALSE)), ]
  go$Term <- factor(go$Term, levels = go$Term)
  
  # 获取每组 top 条目
  top_go <- go %>%
    group_by(ONTOLOGY) %>%
    slice_head(n = topNumber) %>%
    ungroup()
  
  # 构建输出文件路径
  if (!is.null(suffix)) {
    # 添加 suffix 时用_连接
    all_file <- file.path(outputDir, paste0("GOAll_", suffix, ".csv"))
    top_file <- file.path(outputDir, paste0("GOTop", topNumber, "_", suffix, ".csv"))
  } else {
    # 无 suffix 时使用默认文件名
    all_file <- file.path(outputDir, "GOAll.csv")
    top_file <- file.path(outputDir, paste0("GOTop", topNumber, ".csv"))
  }
  
  # 安全写入文件
  tryCatch({
    write.csv(go, all_file, row.names = FALSE)
    write.csv(top_go, top_file, row.names = FALSE)
    message("结果已保存至:\n", all_file, "\n", top_file)
  }, error = function(e) {
    warning("文件保存失败: ", e$message)
  })
  
  return(list(go = go, top_go = top_go))
}


#' KEGG 通路富集分析函数
#'
#' 对基因列表进行 KEGG 通路富集分析，支持自动基因 ID 转换，结果自动保存为 CSV 文件。
#'
#' @param geneDf 数据框，必须包含 ENTREZID 列或 SYMBOL 列（自动转换）
#' @param pAdjustMethod p 值校正方法，默认为 "BH"（Benjamini-Hochberg）
#' @param pvalueCutoff p 值阈值，默认为 1（不过滤）
#' @param qvalueCutoff q 值阈值，默认为 1（不过滤）
#' @param topNumber 显示的最大通路条目数，默认为 30
#' @param suffix 输出文件名的后缀标识（自动添加下划线连接），默认为 NULL
#' @param species 物种选择，可选 "human" 或 "mouse"，默认为 "human"
#' @param outputDir 输出文件保存目录，默认为当前工作目录（"."）
#'
#' @return 返回包含两个元素的列表：
#' \itemize{
#'   \item ekegg - 完整 KEGG 富集结果（包含基因符号的可读格式）
#'   \item top_kegg - 前 topNumber 个显著通路的结果
#' }
#'
#' @examples
#' \dontrun{
#' # 人类基因分析（保存到当前目录）
#' KEGGEnrichment(human_geneDf, suffix = "patient_sample")
#'
#' # 小鼠基因分析（保存到指定目录）
#' KEGGEnrichment(mouse_geneDf, 
#'                species = "mouse", 
#'                outputDir = "./KEGG")
#' }
#'
#' @export
#'
#' @importFrom clusterProfiler enrichKEGG bitr setReadable
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
KEGGEnrichment <- function(geneDf, 
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 1, 
                           qvalueCutoff = 1, 
                           topNumber = 30, 
                           suffix = NULL,
                           species = c("human", "mouse"),
                           outputDir = ".") {
  
  # 隐藏包加载时自动显示的启动消息
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(tidyverse)
    # 根据物种加载对应的 OrgDb
    species <- match.arg(species)
    if (species == "mouse") {
      library(org.Mm.eg.db)
      OrgDb <- org.Mm.eg.db
      organism <- "mmu"
    } else {
      library(org.Hs.eg.db)
      OrgDb <- org.Hs.eg.db
      organism <- "hsa"
    }
  })
  
  # 检查输出目录是否存在，不存在则创建
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
    message("创建输出目录: ", outputDir)
  }
  
  # 参数检查与基因 ID 转换
  if (!"ENTREZID" %in% colnames(geneDf)) {
    if ("SYMBOL" %in% colnames(geneDf)) {
      # 自动从 SYMBOL 转换到 ENTREZID
      convert <- bitr(geneDf$SYMBOL,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = OrgDb)
      
      # 合并转换结果（保留原始数据所有行）
      geneDf <- merge(geneDf, convert, by = "SYMBOL", all.x = TRUE)
      message("已自动将 SYMBOL 转换为 ENTREZID，缺失转换的基因将被过滤")
    } else {
      stop("输入数据框必须包含 ENTREZID 列或 SYMBOL 列")
    }
  }
  
  # 进行 KEGG 富集分析
  ekegg <- enrichKEGG(gene = geneDf$ENTREZID,
                      organism = organism,
                      pAdjustMethod = pAdjustMethod,
                      pvalueCutoff = pvalueCutoff,
                      qvalueCutoff = qvalueCutoff)
  
  ekegg@result$Pathway <- paste(ekegg@result$ID, ekegg@result$Description, sep = ': ')
  ekegg2 <- setReadable(ekegg, OrgDb = OrgDb, keyType = "ENTREZID")
  
  # 提取前 topNumber 个通路
  if (nrow(ekegg2@result) >= topNumber) {
    top_kegg <- ekegg2@result[1:topNumber, ]
  } else {
    top_kegg <- if (nrow(ekegg2) > 0) ekegg2@result else data.frame()
  }
  
  # 构建输出文件路径
  if (!is.null(suffix)) {
    # 添加 suffix 时用_连接
    all_file <- file.path(outputDir, paste0("KEGGAll_", suffix, ".csv"))
    top_file <- file.path(outputDir, paste0("KEGGTop", topNumber, "_", suffix, ".csv"))
  } else {
    # 无 suffix 时使用默认文件名
    all_file <- file.path(outputDir, "KEGGAll.csv")
    top_file <- file.path(outputDir, paste0("KEGGTop", topNumber, ".csv"))
  }
  
  # 安全写入文件
  tryCatch({
    write.csv(ekegg2, all_file, row.names = FALSE)
    write.csv(top_kegg, top_file, row.names = FALSE)
    message("结果已保存至:\n", all_file, "\n", top_file)
  }, error = function(e) {
    warning("文件保存失败: ", e$message)
  })
  
  # 返回结果列表
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
