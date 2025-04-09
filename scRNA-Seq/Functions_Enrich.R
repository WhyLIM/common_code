
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


#' 绘制 GO 和 KEGG 富集分析结果的条形图
#'
#' 该函数接受由 clusterProfiler 生成的 GO/KEGG 富集分析结果数据框
#' 生成带有分类标签的条形图，展示显著富集的通路
#' 按类别分组，支持自定义颜色、排序方式和筛选标准
#'
#' @param go_df GO 富集分析结果数据框，必须包含 ONTOLOGY、Description、p.adjust 等列。可以为 NULL
#' @param kegg_df KEGG 富集分析结果数据框，必须包含 Description、p.adjust 等列。可以为 NULL
#' @param top_n_per_category 每个类别中显示的通路数量。默认为 5
#' @param filter_criterion 用于筛选通路的统计指标，如 "p.adjust"、"pvalue" 或 "qvalue"。默认为 "p.adjust"
#' @param x_axis_column 用于 x 轴的数值列名，可选以下列：
#'   \itemize{
#'     \item `"p.adjust"`：校正后的p值（默认）
#'     \item `"pvalue"`：原始p值
#'     \item `"qvalue"`：校正后的 q 值
#'     \item `"FoldEnrichment"`：富集倍数
#'     \item `"zScore"`：富集 z 分数
#'     \item `"GeneRatio"`：基因比例（如"5/100"会自动计算为 0.05）
#'     \item `"Count"`：差异基因数量
#'   }
#'   数据框中必须包含指定的列。
#' @param x_axis_trans 数据转换方式，可选：
#'   \itemize{
#'     \item `"auto"`：自动选择（默认），对 p 值类（p.adjust/pvalue/qvalue）自动应用 -log10 转换
#'     \item `"-log10"`：强制使用 -log10 转换
#'     \item `"identity"`：不使用转换（原始值）
#'     \item `"log2"`：log2 转换（需正数）
#'   }
#' @param x_axis_label x 轴标签文本。默认为 NULL 时根据`x_axis_column`自动生成：
#'   \itemize{
#'     \item `p.adjust` → "-log10(Adjusted P Value)"
#'     \item `GeneRatio` → "Gene Ratio"
#'     \item 其他列名 → 直接使用列名
#'   }
#'   指定该参数将覆盖自动生成的标签。
#' @param sort_by 排序依据的列名。默认为 "p.adjust"
#' @param filename 输出图片的文件名前缀（不包含扩展名）。默认为 "pathway_plot"
#' @param width 图片宽度（英寸）。默认为 NULL（自动计算）
#' @param height 图片高度（英寸）。默认为 NULL（自动计算）
#' @param custom_colors 自定义颜色向量。默认为 NULL（使用内置配色）
#' @param palette 内置配色方案编号（1-4）。默认为 1
#' @param alpha 颜色透明度（0-1）。默认为 0.8
#' @param bar_width 条形宽度（0-1）。默认为 0.6
#' @param text_size 通路名称文本大小。默认为 5
#' @param gene_text_size 基因 ID 文本大小。默认为 3.5
#' @param save_plot 是否保存图片。默认为 TRUE
#'
#' @return 返回包含两个元素的列表：
#' \itemize{
#'   \item plot - ggplot2 对象，可进一步修改
#'   \item pathway_data - 用于绘图的数据框
#' }
#'
#' @examples
#' \dontrun{
#' # 示例 1：基本用法
#' result <- plot_pathway(go_df = go_results, kegg_df = kegg_results)
#' 
#' # 示例 2：只绘制 KEGG 的结果，使用 GeneRatio 作为 x 轴，显示前 10 条通路
#' result <- plot_pathway(
#'   kegg_df = kegg_results, 
#'   x_axis_column = "GeneRatio",
#'   top_n_per_category = 10,
#'   palette = 2,
#'   alpha = 0.6
#' )
#' 
#' # 查看绘图
#' print(result$plot)
#' 
#' # 查看绘图数据
#' head(result$pathway_data)
#' }
#'
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom gground geom_round_col geom_round_rect
#' @importFrom tibble rowid_to_column
#' @importFrom scales expansion
#' @importFrom grDevices png pdf
plot_pathway <- function(
    go_df = NULL, 
    kegg_df = NULL, 
    top_n_per_category = 5,
    filter_criterion = "p.adjust",
    x_axis_column = "p.adjust",
    x_axis_trans = "auto",
    x_axis_label = NULL,
    sort_by = "p.adjust",
    filename = "pathway_plot",
    width = NULL, 
    height = NULL,
    custom_colors = NULL,
    palette = 1,
    alpha = 0.8,
    bar_width = 0.6,
    text_size = 5,
    gene_text_size = 3.5,
    save_plot = TRUE
) {
  library(dplyr)
  library(ggplot2)
  
  # 检查是否至少提供了一种数据框
  if (is.null(go_df) && is.null(kegg_df)) {
    stop("至少需要提供 GO 或 KEGG 富集分析结果数据框中的一个")
  }
  
  # 检查gground包是否安装
  if (!requireNamespace("gground", quietly = TRUE)) {
    stop("未检测到 gground 包，请使用 devtools::install_github(\"dxsbiocc/gground\") 安装")
  }
  library(gground)
  
  # 获取排序方向
  sort_direction <- if(ascending) 1 else -1
  
  # 构建完整通路数据框
  use_pathway <- NULL
  
  # 需要的列名
  required_cols <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "FoldEnrichment", 
                     "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
  # 处理 GO 数据框
  if (!is.null(go_df)) {
    # 按 ONTOLOGY 分组，筛选每个类别中的通路
    go_filtered <- go_df %>% 
      group_by(ONTOLOGY) %>%
      slice_min(order_by = !!sym(filter_criterion), n = top_n_per_category) %>% 
      ungroup() %>% 
      select(all_of(required_cols))
    # 添加到总数据框
    use_pathway <- go_filtered
  }
  
  # 处理 KEGG 数据框
  if (!is.null(kegg_df)) {
    # 筛选 KEGG 中的通路
    kegg_filtered <- kegg_df %>%
      mutate(ONTOLOGY = 'KEGG') %>%  # 添加 KEGG 类别标记
      slice_min(order_by = !!sym(filter_criterion), n = top_n_per_category) %>% 
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
  category_order <- rev(all_categories)  # 反转顺序
  
  # 处理 pathway 数据框
  use_pathway <- use_pathway %>%
    # 将 ONTOLOGY 转换为因子并指定水平顺序
    mutate(ONTOLOGY = factor(ONTOLOGY, levels = category_order)) %>%
    # 按 ONTOLOGY 和筛选标准排序，为使绘图按 p 降序，此处要升序
    arrange(ONTOLOGY, desc(!!sym(sort_by))) %>%
    # 将 Description 转换为因子，保持当前顺序
    mutate(Description = factor(Description, levels = rev(Description))) %>%
    # 添加行号作为索引
    tibble::rowid_to_column('index')
  
  # 总通路数量，用于确定图形高度
  total_pathways <- nrow(use_pathway)
  
  # 处理 x 轴参数
  if (is.null(x_axis_label)) {
    x_axis_label <- switch(x_axis_column,
                           "p.adjust" = "-log10(Adjusted P Value)",
                           "pvalue" = "-log10(P Value)",
                           "qvalue" = "-log10(Q Value)",
                           "FoldEnrichment" = "Fold Enrichment",
                           "zScore" = "Z Score",
                           "GeneRatio" = "Gene Ratio",
                           "Count" = "Gene Count",
                           x_axis_column  # 默认直接使用列名
    )
  }
  
  # 自动确定数据转换方式
  if (x_axis_trans == "auto") {
    x_axis_trans <- if (x_axis_column %in% c("p.adjust", "pvalue", "qvalue")) "-log10" else "identity"
  }
  
  # 准备 x 轴数据
  if (x_axis_column == "GeneRatio") {
    use_pathway <- use_pathway %>%
      mutate(x_value = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))
  } else {
    use_pathway <- use_pathway %>%
      mutate(x_value = !!sym(x_axis_column))
  }
  
  # 应用数据转换
  if (x_axis_trans == "-log10") {
    use_pathway <- use_pathway %>%
      mutate(x_value = -log10(x_value))
    # 计算 x 轴最大值并添加缓冲
    xaxis_max <- max(use_pathway$x_value, na.rm = TRUE) * 1.1
  } else {
    xaxis_max <- max(use_pathway$x_value, na.rm = TRUE) * 1.1  # 默认增加 10% 空间
  }
  
  # 动态计算 width_factor 以计算左侧分类标签和基因数量点图的宽度
  width_factor <- 1/20 * max(use_pathway$x_value)
  
  # 计算文本位置
  text_positions <- list(
    desc_pos = 1/70 * max(use_pathway$x_value),  # 通路名称位置
    gene_pos = 1/50 * max(use_pathway$x_value)   # 基因 ID 位置
  )
  
  # 准备左侧分类标签数据
  rect.data <- use_pathway %>%
    group_by(ONTOLOGY) %>%
    summarize(n = n(), .groups = "drop") %>%
    ungroup() %>%
    mutate(
      xmin = -2.8 * width_factor,  # 矩形左侧位置
      xmax = -2.2 * width_factor,  # 矩形右侧位置
      ymax = cumsum(n),          # 矩形上边界（累计和）
      ymin = lag(ymax, default = 0) + 0.6,  # 矩形下边界（前一个 ymax 加 0.6）
      ymax = ymax + 0.4          # 调整上边界加 0.4
    )
  
  # 颜色方案
  if (is.null(custom_colors)) {
    # 内置颜色（共 16 种）
    default_colors <- c(
      # 色系 1（红系）
      "#ED6F6E", "#E08D8B", "#F1AEA7", "#EDADC5",
      # 色系 2（橙黄系）
      "#F28147", "#ECB884", "#FAC074", "#FCBB44",
      # 色系 3（蓝绿系）
      "#9FD4AE", "#AAD7C8", "#6CBEC3", "#9ECAE1",
      # 色系 4（蓝紫系）
      "#619CD9", "#92A5D1", "#9D9ECD", "#9584C1"
    )
    
    # 检查 palette 是否有效
    if (!palette %in% 1:4) {
      stop("`palette` 必须为 1 ~ 4 之间的整数")
    }
    
    # 从每个色系中选取第 palette 个颜色
    selected_indices <- seq(palette, by = 4, length.out = 4)
    colors_to_use <- default_colors[selected_indices]
    
    # 根据实际类别数截取颜色
    n_categories <- length(all_categories)
    colors_to_use <- colors_to_use[1:n_categories]
  } else {
    colors_to_use <- custom_colors  # 用户自定义颜色优先
  }
  
  # 创建 ggplot 对象
  p <- ggplot(use_pathway, aes(x = x_value, y = index, fill = ONTOLOGY)) +
    # 添加条形图（圆角矩形）
    gground::geom_round_col(aes(y = Description), width = bar_width, alpha = alpha) +
    scale_y_discrete(limits = rev) +
    # 添加左侧分类标签（圆角矩形）
    gground::geom_round_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
      data = rect.data, radius = unit(2, 'mm'), inherit.aes = FALSE, alpha = alpha
    ) +
    # 添加通路名称文本（左侧对齐）
    geom_text(aes(x = text_positions$desc_pos, label = Description), hjust = 0, size = text_size) +
    # 添加基因 ID 文本（斜体，小字号）
    geom_text(
      aes(x = text_positions$gene_pos, label = geneID, colour = ONTOLOGY), 
      hjust = 0, vjust = 3, size = gene_text_size, fontface = 'italic', 
      show.legend = FALSE
    ) +
    # 添加基因数量点图
    geom_point(aes(x = -width_factor, size = Count), shape = 21, alpha = alpha) +
    # 添加基因数量文本
    geom_text(aes(x = -width_factor, label = Count)) +
    # 设置点的大小范围
    scale_size_continuous(name = 'Count', range = c(5, 16)) +
    # 添加分类标签文本
    geom_text(
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
      data = rect.data, inherit.aes = FALSE, angle = 90
    ) +
    # 设置 y 轴标签为空
    labs(y = NULL) +
    # 设置填充颜色和图例标题
    scale_fill_manual(name = 'Category', values = colors_to_use) +
    # 设置文本颜色
    scale_colour_manual(values = colors_to_use) +
    # 设置 x 轴刻度和扩展
    (if (x_axis_trans == "-log10") {
      scale_x_continuous(
        name = x_axis_label,
        expand = expansion(c(0, 0.05))
      )
    } else {
      scale_x_continuous(
        name = x_axis_label,
        expand = expansion(c(0, 0.05))
      )
    }) +
    # 添加 x 轴线
    annotate("segment", x = 0, y = 0, xend = xaxis_max, yend = 0, linewidth = 1.5) +
    theme_classic() +
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
    ggsave(paste0(filename, ".png"), p, width = width, height = height)
    ggsave(paste0(filename, ".pdf"), p, width = width, height = height)
  }
  
  return(list(plot = p, pathway_data = use_pathway))
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
