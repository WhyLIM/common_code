
library(Seurat)
library(SeuratData)

# Sys.setenv(http_proxy = "http://127.0.0.1:7892")
# Sys.setenv(https_proxy = "http://127.0.0.1:7892")
# InstallData("ifnb")

ifnb <- LoadData("ifnb")
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)

ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))

# ===============================================================================================
setwd("~/Documents/temp")

ego <- readRDS("~/Projects/scLC/results/denovo/W_Penis/me/check-markers-0.5/DEA/Enrich/go_enrichment_results.rds")
ekegg <- readRDS("~/Projects/scLC/results/denovo/W_Penis/me/check-markers-0.5/DEA/Enrich/kegg_enrichment_results.rds")

go_df <- ego$Adipocytes_up$top_go
kegg_df <- ekegg$Adipocytes_up$top_kegg

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
  library(dplyr)
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

plot_pathway(kegg_df, top_n_per_category = 10)


