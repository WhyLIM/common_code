
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


#' 通路富集分析可视化函数
#' 
#' @param go_df GO富集分析结果数据框，可为 NULL
#' @param kegg_df KEGG富集分析结果数据框，可为 NULL
#' @param top_n_per_category 每个类别保留的通路数量，默认为 5
#' @param filter_criterion 用于筛选通路的标准，默认是 "p.adjust"
#' @param second_criterion 相同主要筛选标准时的次要筛选标准，默认是 "Count"
#' @param filename 输出文件名，默认为 "pathway_plot"
#' @param width 输出图形宽度，默认为 NULL（自动计算）
#' @param height 输出图形高度，默认为 NULL（自动计算）
#' @param custom_colors 自定义颜色方案，默认为 NULL（使用内置颜色）
#' @param sort_by 排序依据，默认为 "p.adjust"
#' @param ascending 是否升序排列，默认为 TRUE（p.adjust 值小的在前）
#' @param bar_width 条形图宽度，默认为 0.6
#' @param text_size 通路描述文本大小，默认为 5
#' @param gene_text_size 基因 ID 文本大小，默认为 3.5
#' @param save_plot 是否保存绘图文件，默认为 TRUE
#' @return ggplot对象
#' 
#' @examples
#' # 基本用法 - 使用默认参数
#' plot_pathway(go_df, kegg_df)
#' 
#' 只展示 GO 通路，不展示 KEGG 通路
#' plot_pathway(go_df = go_df, kegg_df = NULL)
#' 
#' 每个类别选择 10 个通路
#' plot_pathway(go_df, kegg_df, top_n_per_category = 10)
#' 
#' 使用自定义文件名和尺寸
#' plot_pathway(go_df, kegg_df, width = 12, height = 10)
#' 
#' 使用其他筛选标准（如使用 pvalue 而不是 p.adjust）
#' plot_pathway(go_df, kegg_df, filter_criterion = "pvalue", sort_by = "pvalue")
plot_pathway <- function(
    go_df = NULL, 
    kegg_df = NULL, 
    top_n_per_category = 5,
    filter_criterion = "p.adjust",
    second_criterion = "Count",
    filename = "pathway_plot",
    width = NULL, 
    height = NULL,
    custom_colors = NULL,
    sort_by = "p.adjust",
    ascending = TRUE,
    bar_width = 0.6,
    text_size = 5,
    gene_text_size = 3.5,
    save_plot = TRUE
) {
  # 检查是否至少提供了一种数据框
  if (is.null(go_df) && is.null(kegg_df)) {
    stop("至少需要提供 GO 或 KEGG 富集分析结果数据框中的一个")
  }
  
  # 获取排序方向
  sort_direction <- if(ascending) 1 else -1
  
  # 构建完整通路数据框
  use_pathway <- NULL
  
  # 需要的列名
  required_cols <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", 
                     "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  # 处理 GO 数据框
  if (!is.null(go_df)) {
    # 按 ONTOLOGY 分组，选取每个类别中筛选标准值最佳的通路
    go_filtered <- go_df %>% 
      group_by(ONTOLOGY) %>%
      # 使用动态筛选
      top_n(top_n_per_category, wt = sort_direction * !!sym(filter_criterion)) %>%
      # 二级筛选标准
      group_by(!!sym(filter_criterion)) %>%
      top_n(1, wt = !!sym(second_criterion)) %>%
      ungroup() %>% 
      select(all_of(required_cols))
    # 添加到总数据框
    use_pathway <- go_filtered
  }
  
  # 处理 KEGG 数据框
  if (!is.null(kegg_df)) {
    # 选取 KEGG 中筛选标准值最佳的通路
    kegg_filtered <- kegg_df %>%
      top_n(top_n_per_category, wt = sort_direction * !!sym(filter_criterion)) %>%
      # 二级筛选标准
      group_by(!!sym(filter_criterion)) %>%
      top_n(1, wt = !!sym(second_criterion)) %>%
      mutate(ONTOLOGY = 'KEGG') %>%  # 添加 KEGG 类别标记
      ungroup() %>% 
      select(all_of(required_cols))
    # 添加到总数据框
    if (is.null(use_pathway)) {
      use_pathway <- kegg_filtered
    } else {
      use_pathway <- rbind(use_pathway, kegg_filtered)
    }
  }
  
  # 获取所有类别
  all_categories <- unique(use_pathway$ONTOLOGY)
  category_order <- rev(all_categories)  # 默认反转顺序，使 GO 类别在上，KEGG 在下
  
  # 处理 pathway 数据框
  use_pathway <- use_pathway %>%
    # 将 ONTOLOGY 转换为因子并指定水平顺序
    mutate(ONTOLOGY = factor(ONTOLOGY, levels = category_order)) %>%
    # 按 ONTOLOGY 和筛选标准排序
    arrange(ONTOLOGY, !!sym(sort_by)) %>%
    # 将 Description 转换为因子，保持当前顺序
    mutate(Description = factor(Description, levels = Description)) %>%
    # 添加行号作为索引
    tibble::rowid_to_column('index')
  
  # 总通路数量，用于确定图形高度
  total_pathways <- nrow(use_pathway)
  
  # 左侧分类标签和基因数量点图的宽度
  width_factor <- 0.5
  
  # 计算 x 轴最大值（基于 -log10(p.adjust) 的最大值加 1）
  if ("p.adjust" %in% colnames(use_pathway)) {
    xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
  } else {
    # 如果没有 p.adjust 列，使用默认值
    xaxis_max <- 10
  }
  
  # 准备左侧分类标签数据
  rect.data <- use_pathway %>%
    group_by(ONTOLOGY) %>%
    summarize(n = n(), .groups = "drop") %>%
    ungroup() %>%
    mutate(
      xmin = -3 * width_factor,  # 矩形左侧位置
      xmax = -2 * width_factor,  # 矩形右侧位置
      ymax = cumsum(n),          # 矩形上边界（累计和）
      ymin = lag(ymax, default = 0) + 0.6,  # 矩形下边界（前一个 ymax 加 0.6）
      ymax = ymax + 0.4          # 调整上边界加 0.4
    )
  
  # 自定义颜色方案
  if (is.null(custom_colors)) {
    default_colors <- c(
      "#9584C1", "#FCBB44", "#ED6F6E", "#6CBEC3", 
      "#FED98E", "#F28147", "#619CD9", "#EDADC5", 
      "#9ECAE1", "#ECB884", "#92A5D1", "#E08D8B", 
      "#AAD7C8", "#F1AEA7", "#AF8CBB", "#FAC074", 
      "#9D9ECD", "#9FD4AE", "#D9B9D4"
    )
    # 确保颜色数量足够
    colors_to_use <- default_colors[1:length(all_categories)]
  } else {
    colors_to_use <- custom_colors
  }
  
  # 检查是否加载了 gground 包（用于 geom_round_col 和 geom_round_rect）
  if (!requireNamespace("gground", quietly = TRUE)) {
    message("未检测到 gground 包，尝试使用 geom_col 和 geom_rect 替代...")
    # 创建 ggplot 对象（无圆角）
    p <- ggplot(use_pathway, aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
      # 添加条形图（使用普通矩形）
      geom_col(aes(y = Description), width = bar_width, alpha = 0.8)
    
    # 添加左侧分类标签（矩形，无圆角）
    p <- p + geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
      data = rect.data, inherit.aes = FALSE
    )
  } else {
    # 创建 ggplot 对象（带圆角）
    p <- ggplot(use_pathway, aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
      # 添加条形图（使用圆角矩形）
      gground::geom_round_col(aes(y = Description), width = bar_width, alpha = 0.8)
    
    # 添加左侧分类标签（圆角矩形）
    p <- p + gground::geom_round_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
      data = rect.data, radius = unit(2, 'mm'), inherit.aes = FALSE
    )
  }
  
  # 添加其他图层
  p <- p +
    # 添加通路名称文本（左侧对齐）
    geom_text(aes(x = 0.05, label = Description), hjust = 0, size = text_size) +
    # 添加基因 ID 文本（斜体，小字号）
    geom_text(
      aes(x = 0.1, label = geneID, colour = ONTOLOGY), 
      hjust = 0, vjust = 3, size = gene_text_size, fontface = 'italic', 
      show.legend = FALSE
    ) +
    # 添加基因数量点图
    geom_point(aes(x = -width_factor, size = Count), shape = 21) +
    # 添加基因数量文本
    geom_text(aes(x = -width_factor, label = Count)) +
    # 设置点的大小范围
    scale_size_continuous(name = 'Count', range = c(5, 16)) +
    # 添加分类标签文本
    geom_text(
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
      data = rect.data, inherit.aes = FALSE
    ) +
    # 添加 x 轴线
    geom_segment(
      aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
      linewidth = 1.5, inherit.aes = FALSE
    ) +
    # 设置 y 轴标签为空
    labs(y = NULL) +
    # 设置填充颜色和图例标题
    scale_fill_manual(name = 'Category', values = colors_to_use) +
    # 设置文本颜色
    scale_colour_manual(values = colors_to_use) +
    # 设置x轴刻度和扩展
    scale_x_continuous(
      name = "-log10(p.adjust)",
      breaks = seq(0, ceiling(xaxis_max), 2), 
      expand = expansion(c(0, 0))
    ) +
    # 使用经典主题
    theme_classic() +
    # 自定义主题设置
    theme(
      axis.text.y = element_blank(),  # 隐藏 y 轴文本
      axis.line = element_blank(),    # 隐藏轴线
      axis.ticks.y = element_blank(), # 隐藏 y 轴刻度
      legend.title = element_text()   # 设置图例标题文本
    )
  
  # 自动计算合适的图形尺寸
  if (is.null(width)) {
    width <- 10  # 默认宽度
  }
  if (is.null(height)) {
    # 根据通路数量动态计算高度
    height <- total_pathways * 0.55  # 每个通路 0.55 英寸高度
  }
  
  if (save_plot) {
    # 保存图形
    ggsave(paste0(filename, ".png"), p, width = width, height = height)
    ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
  }
  
  return(p)
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
