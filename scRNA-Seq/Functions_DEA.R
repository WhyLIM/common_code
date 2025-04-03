
############################################ 执行 PCA 的函数 #############################################
pcaPlot <- function(limma_exp, group, palette = NULL) {
  library(FactoMineR)
  library(factoextra)
  library(tidyverse)
  # 画 PCA 图时要求是行名为样本名，列名为基因。转置表达矩阵
  pca_exp <- t(limma_exp)
  
  # 将表达矩阵转换为数据框，并合并分组信息
  pca_exp_with_group <- pca_exp %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "id") %>%   # 将行名转换为列
    merge(., group, by = "id") %>%       # 合并 group 信息
    column_to_rownames(var = "id")       # 将 id 列转换回行名
  
  # 进行PCA分析
  pca_data <- PCA(pca_exp, graph = FALSE)
  
  # 绘制PCA图
  fviz_pca_ind(pca_data,
               geom.ind = "point",
               col.ind = pca_exp_with_group$group, # 根据分组上色
               palette = palette,  # 用户可以自定义调色板
               addEllipses = TRUE, # 添加浓度椭圆
               legend.title = "Groups"
  )
}


######################################## 带有去批次的差异分析函数 ########################################
# group 参数必须含有名为 group 的列
runLimma <- function(exp, group, normalize = TRUE, remove_batch_effect = FALSE, batch = NULL, batch_design = NULL, directional = FALSE) {
  library(limma)
  library(edgeR)
  library(tidyverse)
  
  # 分组信息变为因子向量
  conditions <- factor(group$group)
  
  # 设计矩阵
  design <- model.matrix(~0 + conditions)
  # design <- model.matrix(~conditions)
  colnames(design) <- levels(conditions)
  print(design)
  
  if (normalize) {
    # RNA-Seq 数据归一化，DGEList 是为整数 count 矩阵设计的
    dge <- DGEList(counts = as.matrix(exp), group = conditions)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = TRUE)
    exprs <- v$E
  } else {
    exprs <- as.matrix(exp)
  }
  
  # 去除批次效应
  pca_plot <- NA
  if (remove_batch_effect) {
    if (is.null(batch)) {
      stop("Batch information must be provided if remove_batch_effect is TRUE")
    }
    exprs <- removeBatchEffect(exprs, batch = batch, design = batch_design)
    pca_plot <- pcaPlot(exprs, group, palette = c("#00AFBB", "#E7B800"))
  }
  
  # 线性模型和经验贝叶斯
  fit <- lmFit(exprs, design)
  fit <- eBayes(fit)
  
  # 自动生成所有可能的对比矩阵
  unique_conditions <- levels(conditions)
  num_conditions <- length(unique_conditions)
  contrast_str <- c()
  for (i in 1:(num_conditions - 1)) {
    for (j in (i + 1):num_conditions) {
      contrast_name <- paste0(unique_conditions[i], "-", unique_conditions[j])
      contrast_str <- c(contrast_str, contrast_name)
      if (directional) {
        contrast_name_reverse <- paste0(unique_conditions[j], "-", unique_conditions[i])
        contrast_str <- c(contrast_str, contrast_name_reverse)
      }
    }
  }
  contrasts <- makeContrasts(contrasts = contrast_str, levels = colnames(design))
  
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  
  # 生成详细的结果总结
  summary_results <- list()
  for (contrast_name in contrast_str) {
    tT <- topTable(fit2, coef = contrast_name, adjust = "fdr", sort.by = "B", number = nrow(exprs))
    regulated <- ifelse(tT$logFC > 0, "Up", "Down")
    tT$Regulation <- regulated
    lfcs <- c(log2(1.2), log2(1.3), log2(1.5), 1)
    all.deg.cnt <- cbind()
    for (lfc in lfcs) {
      deg1 <- regulated[which(abs(tT$logFC) > lfc & tT$P.Value < 0.05)]
      deg2 <- regulated[which(abs(tT$logFC) > lfc & tT$P.Value < 0.01)]
      deg3 <- regulated[which(abs(tT$logFC) > lfc & tT$adj.P.Val < 0.1)]
      deg4 <- regulated[which(abs(tT$logFC) > lfc & tT$adj.P.Val < 0.05)]
      deg5 <- regulated[which(abs(tT$logFC) > lfc & tT$adj.P.Val < 0.01)]
      all.deg.cnt <- cbind(all.deg.cnt, c(
        paste0(sum(deg1 == "Up"), "|", sum(deg1 == "Down")),
        paste0(sum(deg2 == "Up"), "|", sum(deg2 == "Down")),
        paste0(sum(deg3 == "Up"), "|", sum(deg3 == "Down")),
        paste0(sum(deg4 == "Up"), "|", sum(deg4 == "Down")),
        paste0(sum(deg5 == "Up"), "|", sum(deg5 == "Down"))
      ))
    }
    row.names(all.deg.cnt) <- c("p<0.05", "p<0.01", "FDR<0.1", "FDR<0.05", "FDR<0.01")
    colnames(all.deg.cnt) <- paste0(c("1.2", "1.3", "1.5", "2"), "-fold")
    summary_results[[contrast_name]] <- list(DEG = tT, Summary = all.deg.cnt)
  }
  
  return(list(
    CombinedFit = fit2,
    DetailedResults = summary_results, 
    ExpressionMatrix = exprs,
    PCAPlot = pca_plot
  ))
}


############################################# 获得总结的函数 #############################################
getSummary <- function(results, contrast_name) {
  if (!contrast_name %in% names(results$DetailedResults)) {
    stop("The specified contrast is not found in the results.")
  }
  
  summary <- results$DetailedResults[[contrast_name]]$Summary
  print(summary)
}

########################################### 获得差异基因的函数 ###########################################
getDEGs <- function(results, contrast_name, pvalue_threshold = 0.05, fdr_threshold = 0.05, logFC_threshold = 1) {
  if (!contrast_name %in% names(results$DetailedResults)) {
    stop("The specified contrast is not found in the results.")
  }
  
  tT <- results$DetailedResults[[contrast_name]]$DEG
  if (is.null(pvalue_threshold) && is.null(fdr_threshold) && is.null(logFC_threshold)) {
    degs <- tT
  } else {
    degs <- tT[tT$P.Value < pvalue_threshold & tT$adj.P.Val < fdr_threshold & abs(tT$logFC) > logFC_threshold, ]
  }
  
  return(degs)
}


############################################ 绘制火山图的函数 ############################################
volcanoPlot <- function(data, logFC_col, pvalue_col, title = "", 
                        colors = c("#0072B5", "grey", "#BC3C28"), labels = NULL, 
                        pvalue_threshold = 0.05, logFC_threshold = 1) {
  library(ggrepel)
  library(dplyr)
  library(rlang)
  
  # 检查颜色向量长度
  if (length(colors) != 3) {
    stop("Please provide exactly three colors.")
  }
  
  # 将列名字符串转换为 symbols
  log2FC_sym <- sym(logFC_col)
  pvalue_sym <- sym(pvalue_col)
  
  # 在函数内部计算 Significance 列
  data <- data %>%
    mutate(Significance = ifelse(
      !!pvalue_sym < pvalue_threshold & abs(!!log2FC_sym) > logFC_threshold,
      ifelse(!!log2FC_sym > logFC_threshold, "Up", "Down"),
      "NS"
    ))
  
  # 如果提供了 labels，添加一个新列用于标注
  if (!is.null(labels)) {
    data$Label <- ifelse(rownames(data) %in% labels, 1, 0)
  } else {
    data$Label <- 0
  }
  
  # 定义主题
  mytheme <- theme_bw() + 
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 13), 
          axis.title = element_text(size = 13))
  
  # 绘图
  p <- ggplot(data, aes(x = !!log2FC_sym, y = -log10(!!pvalue_sym), color = Significance)) +
    geom_point(alpha = 0.5) + 
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = colors) +
    labs(title = title, x = "log2(FoldChange)", y = paste0("-log10(", pvalue_col, ")")) +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), 
               colour = "black", linetype = "dashed") +
    geom_hline(yintercept = -log10(pvalue_threshold), 
               colour = "black", linetype = "dashed") +
    mytheme
  
  # 如果有标签，添加标签
  if (!is.null(labels)) {
    p <- p + geom_text_repel(aes(label = ifelse(Label == 1, rownames(data), "")), 
                             arrow = arrow(ends="first", length = unit(0.01, "npc")), show.legend = F)
  }
  
  return(p)
}


########################################### 绘制火山图的函数2 ############################################
volcanoPlot2 <- function(limma_result, logfc_col = "logFC", pvalue_col = "adj.P.Val", 
                         significance_col = "Significance", 
                         colors = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                         title = "", pvalue_threshold = 0.05, logFC_threshold = 1, 
                         labels = NULL, legend_position = "right") {
  
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(grid)
  
  # 检查数据有效性
  if (!all(c(logfc_col, pvalue_col) %in% colnames(limma_result))) {
    stop("logfc_col or pvalue_col not found in dataframe")
  }
  
  # 动态计算显著性列
  if (!significance_col %in% colnames(limma_result)) {
    limma_result <- limma_result %>%
      mutate(!!sym(significance_col) := case_when(
        !!sym(pvalue_col) < pvalue_threshold & !!sym(logfc_col) > logFC_threshold ~ "Up",
        !!sym(pvalue_col) < pvalue_threshold & !!sym(logfc_col) < -logFC_threshold ~ "Down",
        TRUE ~ "NS"
      ))
  }
  
  # 准备标注逻辑
  if (!is.null(labels)) {
    valid_labels <- intersect(labels, rownames(limma_result))
    if (length(valid_labels) == 0) warning("No valid labels found")
    limma_result$custom_label <- ifelse(rownames(limma_result) %in% valid_labels, rownames(limma_result), "")
  }
  
  # 核心绘图
  p <- ggplot(limma_result, aes(x = !!sym(logfc_col), 
                                y = -log10(!!sym(pvalue_col)),
                                color = !!sym(logfc_col))) +
    
    # 散点图层
    geom_point(aes(size = -log10(!!sym(pvalue_col))), alpha = 0.8) +
    
    # 颜色梯度映射
    scale_color_gradientn(
      colours = colors,
      values = scales::rescale(seq(-max(abs(limma_result[[logfc_col]])), 
                                   max(abs(limma_result[[logfc_col]])), 
                                   length.out = 5)),
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5)
    ) +
    
    # 自定义标注
    if (!is.null(labels)) {
      ggrepel::geom_text_repel(
        aes(label = custom_label),
        color = "black",
        size = 4,
        fontface = "bold"
      )
    }
  
  # 添加图形辅助元素
  p <- p +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dotted") +
    scale_size_continuous(range = c(1, 5)) +
    labs(
      title = title,
      x = expression(log[2]("Fold Change")),
      y = bquote(-log[10](.(pvalue_col))),
      color = "Log2 Fold Change"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = legend_position,
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # 添加自定义箭头标注
  p <- p +
    annotation_custom(
      grob = segmentsGrob(
        x0 = unit(0.1, "npc"), x1 = unit(0.3, "npc"),
        y0 = unit(0.95, "npc"), y1 = unit(0.95, "npc"),
        arrow = arrow(type = "closed", length = unit(2, "mm")),
        gp = gpar(col = "#00cec9", lwd = 2)
      )
    ) +
    annotation_custom(
      grob = textGrob(
        label = "Down-regulated",
        x = unit(0.2, "npc"), 
        y = unit(0.92, "npc"),
        gp = gpar(col = "#00cec9", fontface = 2)
      )
    ) +
    annotation_custom(
      grob = segmentsGrob(
        x0 = unit(0.9, "npc"), x1 = unit(0.7, "npc"),
        y0 = unit(0.95, "npc"), y1 = unit(0.95, "npc"),
        arrow = arrow(type = "closed", length = unit(2, "mm")),
        gp = gpar(col = "#fdcb6e", lwd = 2)
      )
    ) +
    annotation_custom(
      grob = textGrob(
        label = "Up-regulated",
        x = unit(0.8, "npc"),
        y = unit(0.92, "npc"),
        gp = gpar(col = "#fdcb6e", fontface = 2)
      )
    )
  
  print(p)
}