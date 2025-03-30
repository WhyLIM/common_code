
############################################## 单细胞分析函数 ############################################
# 初步过滤的函数
basic_qc <- function(seurat_obj, dir = NULL) {
  # 计算特定基因组的比例（线粒体、核糖体、血红蛋白）
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^mt-", col.name = "percent_mito") %>% 
    PercentageFeatureSet(., pattern = "^Rp[sl]", col.name = "percent_ribo") %>% 
    PercentageFeatureSet(., pattern = "^Hb[^p]", col.name = "percent_hb")
  
  # 可视化线粒体、核糖体、血红蛋白基因比例
  p <- VlnPlot(seurat_obj, group.by = "orig.ident", 
               features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb"), 
               pt.size = 0, ncol = 3)
  ggsave(file.path(dir, "Vlnplot_mito_ribo_hb.png"), p, width = length(unique(seurat_obj$orig.ident))/2 + 10, height = 10)
  ggsave(file.path(dir, "Vlnplot_mito_ribo_hb.pdf"), p, width = length(unique(seurat_obj$orig.ident))/2 + 10, height = 10)
  
  return(seurat_obj)
}
# 用法示例：
# seurat_object <- basic_qc(seurat_object)


# 绘制组间表达差异配对箱线图的函数
plot_boxpair <- function(seurat_object, markers, target_genes, group = "group", save_plot = TRUE, 
                         filename = "boxpair", title = "", width = 4, height = 6, group_colors = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(ggrepel)
  
  # 1. 筛选在 seurat_object 中存在的目标基因
  valid_genes <- intersect(target_genes, rownames(seurat_object))
  if (length(valid_genes) == 0) {
    stop("在 seurat 对象中未找到任何目标基因！")
  }
  
  # 2. 筛选出在 markers 中差异显著（p_val_adj < 0.05）的基因
  degs <- markers %>% filter(p_val_adj < 0.05)
  target_degs <- intersect(valid_genes, rownames(degs))
  
  # 3. 提取表达矩阵及相关元数据（只使用 group 信息）
  expr_data <- FetchData(seurat_object, vars = c(valid_genes, group))
  
  # 4. 计算每个组内各目标基因的平均表达值，并转换为长格式数据
  avg_exp <- expr_data %>%
    group_by(across(all_of(group))) %>%
    summarise(across(all_of(valid_genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(cols = all_of(valid_genes),
                 names_to = "gene",
                 values_to = "expression")
  
  # 5. 绘制成对图
  p <- ggpaired(avg_exp, x = group, y = "expression", color = group, palette = "npg", point.size = 3) +
    stat_compare_means(aes_string(x = group, y = "expression"), paired = TRUE, method = "t.test", size = 4) +
    geom_text_repel(aes(label = gene), size = 4, box.padding = 0.5, point.padding = 0.3, segment.color = "grey50",
                    color = ifelse(avg_exp$gene %in% target_degs, "#8e44ad", "black")) +
    labs(title = title, x = "Group", y = "Average Expression", color = "Group") +
    theme_classic()
  
  if (!is.null(group_colors)) {
    p <- p + scale_color_manual(values = group_colors)
  }
  
  # 显示图形
  print(p)
  
  # 6. 如果需要保存图形，则输出为 PDF 和 PNG 文件
  if (save_plot) {
    ggsave(filename = paste0(filename, ".pdf"), plot = p, width = width, height = height)
    ggsave(filename = paste0(filename, ".png"), plot = p, width = width, height = height)
  }
  
}
# 用法示例：
# group_colors <- c(D = "#dd5633", L = "#5fb9d4")
# plot_boxpair(seurat_object, markers = markers, target_genes = genes_CR,
#              group = "group", save_plot = T, filename = "boxpair_CR", 
#              width = 6, height = 8, group_colors = group_colors)


#' boxpair_gene_set
#'
#' 绘制指定基因集在多个细胞类型中的组间表达差异配对箱线图，并将图像保存为 PDF 和 PNG 文件。
#'
#' @param seurat_list      一个包含多个 Seurat 对象的列表，每个代表一个细胞类型。
#' @param deg_list         一个与 seurat_list 对应的 DEG（差异表达）数据框列表，行为基因，必须包含 `p_val_adj` 列。
#' @param genes            目标基因组成的字符向量，用于绘图。
#' @param output_dir       字符串，指定图像保存的输出目录（建议为基因集命名的子文件夹，如 "boxpair_plots/CR"）。
#' @param group            Seurat metadata 中的分组变量（例如 `"condition"` 或 `"group"`）。
#' @param group_colors     一个命名颜色向量，例如 `c("Control" = "#4DBBD5", "Treatment" = "#E64B35")`，用于手动指定分组颜色。
#' @param base_height      数值型，控制图像高度的基准值（单位：英寸），每个基因会增加一定高度。
#' @param base_width       数值型，图像的宽度（单位：英寸）。
#' @param save_plot        逻辑值，是否保存图像文件，默认为 TRUE。
#'
#' @return 无返回值（invisible NULL）。函数的主要作用是保存绘图文件。
#'
#' @examples
#' gene_sets <- list(
#'   CR = list(genes = genes_CR, output_dir = "boxpair_plots/CR", group = "mfield", group_colors = group_colors),
#'   FS = list(genes = genes_FS, output_dir = "boxpair_plots/FS", group = "mfield", group_colors = group_colors)
#' )
#'
#' # 批量绘图
#' lapply(gene_sets, function(x) {
#'   boxpair_gene_set(
#'     seurat_list = seurat_cells,
#'     deg_list = deg_list,
#'     genes = x$genes,
#'     output_dir = x$output_dir,
#'     group = x$group,
#'     group_colors = x$group_colors,
#'     base_height = 0.5,
#'     base_width = 6
#'   )
#' })
#'
boxpair_gene_set <- function(seurat_list, deg_list, genes, 
                             output_dir = "boxpair_plots", 
                             group = "group", 
                             group_colors = NULL, 
                             base_height = 0.5, 
                             base_width = 6,
                             save_plot = TRUE) {
  
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(ggrepel)
  
  # 自动提取基因集名称
  gene_set_name <- basename(normalizePath(output_dir, mustWork = FALSE))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (celltype in names(seurat_list)) {
    seurat_obj <- seurat_list[[celltype]]
    markers <- deg_list[[celltype]]
    
    valid_genes <- intersect(genes, rownames(seurat_obj))
    if (length(valid_genes) == 0) {
      message(paste0("[跳过] ", celltype, " 中无有效基因"))
      next
    }
    
    degs <- markers %>% filter(p_val_adj < 0.05)
    target_degs <- intersect(valid_genes, rownames(degs))
    
    expr_data <- FetchData(seurat_obj, vars = c(valid_genes, group))
    
    avg_exp <- expr_data %>%
      group_by(across(all_of(group))) %>%
      summarise(across(all_of(valid_genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      pivot_longer(cols = all_of(valid_genes), names_to = "gene", values_to = "expression")
    
    file_prefix <- paste0(output_dir, "/boxpair_", gene_set_name, "_", gsub("[/ ]", "_", celltype))
    title <- paste("Expression changes of", gene_set_name, "genes in", celltype)
    height <- base_height * length(valid_genes) + 2
    
    p <- ggpaired(avg_exp, x = group, y = "expression", color = group, 
                  point.size = 3, palette = if (is.null(group_colors)) "npg" else NULL) +
      stat_compare_means(aes_string(x = group, y = "expression"), paired = TRUE, method = "t.test", size = 4) +
      geom_text_repel(aes(label = gene), size = 4, box.padding = 0.5, point.padding = 0.3, segment.color = "grey50",
                      color = ifelse(avg_exp$gene %in% target_degs, "#8e44ad", "black")) +
      labs(title = title, x = "Group", y = "Average Expression", color = "Group") +
      theme_classic()
    
    if (!is.null(group_colors)) {
      p <- p + scale_color_manual(values = group_colors)
    }
    
    if (save_plot) {
      ggsave(filename = paste0(file_prefix, ".pdf"), plot = p, width = base_width, height = height)
      ggsave(filename = paste0(file_prefix, ".png"), plot = p, width = base_width, height = height)
    }
  }
}


# 绘制基因集的基因表达图的函数
featureplot_gene_set <- function(seurat_obj, 
                                 genes, 
                                 output_dir, 
                                 file_prefix, 
                                 base_width = 2, 
                                 base_height = 2, 
                                 ncol = 3) {
  
  # 过滤存在于 Seurat 对象中的基因
  valid_genes <- intersect(genes, rownames(seurat_obj))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (length(valid_genes) == 0) {
    message(paste("No valid genes found for", file_prefix))
    return(NULL)
  }
  
  filename <- paste0(output_dir, "/featureplot_", file_prefix)
  title <- paste("Expression of", file_prefix, "genes")
  
  # 计算行数和默认尺寸
  nrow <- ceiling(length(valid_genes) / ncol)
  
  # 创建绘图对象
  p <- sc_feature(seurat_obj, features = valid_genes, ncol = ncol) +
    labs(title = title)
  
  # 计算宽度和高度
  height <- base_height * nrow
  width <- base_width * ncol
  
  # 保存图片
  ggsave(paste0(filename, ".png"), plot = p, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(filename, ".pdf"), plot = p, width = width, height = height, limitsize = FALSE)
  
}


# 绘制 UMAP/T-SNE 的函数
plot_seurat_dim <- function(seurat_obj, 
                            reduction = "umap", 
                            ident = "orig.ident", 
                            pointsize = 1, 
                            height = 6, 
                            width = 7, 
                            filename = "Dimplot",
                            custom_colors = NULL,  # 长度大于 ident 数量即可
                            name_to_color = NULL,  # 直接为每个 ident 指定颜色
                            label_color = '#5D478B',
                            legend_size = 5,
                            color_palette = 'H',
                            title = "",
                            legend_title = ident, 
                            save_plot = TRUE) { 
  library(ggsc)
  library(ggrepel)
  
  # 设置身份标识符
  Idents(seurat_obj) <- ident
  ids <- Idents(seurat_obj)
  unique_ids <- unique(ids)
  
  # 自动颜色映射逻辑
  if (!is.null(name_to_color)) {
    color_scale <- scale_color_manual(values = name_to_color)
    color_map <- as.list(name_to_color)
  } else if (!is.null(custom_colors)) {
    if (length(unique_ids) > length(custom_colors)) {
      stop("颜色数量少于簇的数量！")
    }
    color_map <- setNames(custom_colors, unique_ids)
    color_scale <- scale_color_manual(values = color_map)
  } else {
    color_map <- viridis::viridis_pal(option = color_palette)(length(unique_ids))
    names(color_map) <- unique_ids
    color_scale <- scale_color_manual(values = color_map)
  }
  
  # 生成可视化图
  p <- sc_dim(seurat_obj, reduction = reduction, pointsize = pointsize) + 
    sc_dim_geom_label(geom = ggrepel::geom_text_repel, 
                      color = label_color, 
                      bg.color = 'white', 
                      nudge_x = 0, nudge_y = 0,      # 可选：微调标签位置
                      size = legend_size - 1         # 可选：调整标签字体大小
    ) + 
    ggtitle(title) + 
    color_scale + 
    guides(color = guide_legend(
      title = legend_title,  
      override.aes = list(size = legend_size)
    ))
  
  # 保存文件
  if (save_plot) {
    ggsave(paste0(filename, ".pdf"), plot = p, height = height, width = width)
    ggsave(paste0(filename, ".png"), plot = p, height = height, width = width)
  }
  
  # 返回颜色映射关系
  return(list(plot = p, color_map = color_map))
}
# 用法示例：
# plot_seurat_dim(seurat_object, ident = "RNA_snn_res.0.8", filename = "UMAP_cluster")
# 优先级顺序：name_to_color > custom_colors >  viridis 调色板
# name_to_color = c("Fibroblasts" = "blue", "Other Cells" = "#FF5733")
# plot_seurat_dim(seurat_object, reduction = "umap", ident = "RNA_snn_res.0.8", 
#                 height = 6, width = 7, filename = "UMAP_cluster",
#                 custom_colors = NULL, name_to_color = NULL, 
#                 label_color = '#5D478B', legend_size = 5, color_palette = 'H',
#                 title = "", legend_title = ident)


# 检测每个 marker 基因集的表达情况（Dotplot）的函数
plot_markers <- function(seurat_obj, marker_list, group.by) {
  library(ggsc)
  
  lapply(names(marker_list), function(cell_type) {
    # 提取该细胞类型所有标记基因
    genes <- marker_list[[cell_type]]
    
    # 验证基因是否存在
    valid_genes <- genes[genes %in% rownames(seurat_obj)]
    if (length(valid_genes) == 0) {
      message("Skipping ", cell_type, ": no valid genes found.")
      return(NULL)
    }
    
    # 绘图
    p <- sc_dot(
      object = seurat_obj,
      features = valid_genes,
      group.by = group.by,
      dot.scale = 6
    ) + 
      scale_color_gradient(low = "#4d8076", high = "red") +
      ggtitle(cell_type) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
      )
    
    # 动态尺寸调整规则
    n_genes <- length(valid_genes)
    base_width <- max(3, n_genes * 1.2)  # 每个基因至少 1.5 英寸宽度
    base_height <- max(4, n_genes * 0.3) # 高度与基因数量相关
    
    # 保存PDF
    ggsave(
      filename = paste0("check_", gsub("[/]", "_", cell_type), ".pdf"),
      plot = p,
      width = base_width,
      height = base_height,
      limitsize = FALSE
    )
    
    return(p)
  })
}
# 使用示例：
# marker_list 的格式：
# marker_genes <- list(
#   # T 细胞
#   T_Cell = c("Cd3g","Cd3e","Cd4","Cd8a","Klrd1"),
#   # NK 细胞
#   NK_Cell = c("Ncr1", "Klrb1a", "Klrk1", "Gzmb", "Prf1"),
#   # 浆细胞
#   Plasma = c('Sdc1', 'Prdm1', 'Xbp1', 'Tnfrsf17'),
#   # B 细胞
#   B_Cell = c("Cd19", "Cd79a", "Cd79b", "Ms4a1"),
#   # 单核/巨噬细胞
#   Mac_Mono = c("Adgre1","Csf1r","Trem2", "Msr1"),
#   # 树突状细胞
#   DCs = c("Siglech", "Clec9a", "Cd209a", "Flt3"),
#   # 中性粒细胞
#   Neu = c("S100a8","S100a9",'Ly6g', 'Cxcr2', 'Cd177'),
#   # 肥大细胞
#   Mast_Cell = c("Fcer1a", "Ms4a2", "Cpa3", "Mcpt8"),
#   # 内皮细胞
#   Endothelial_Cells = c("Kdr","Emcn","Pecam1","Sparc"),
#   # 红系祖细胞
#   Erythroid_Progenitor = c("Gata1","Klf1", "Epor", "Slc4a1", "Tspo2"),
#   # 成熟红细胞
#   Erythrocyte = c("Hba-a1","Hba-a2","Gypa")
# )
# plot_markers(seurat_object, marker_genes, group.by = "RNA_snn_res.0.8")


# 绘制基因集的小提琴图的函数
vlnplot_gene_set <- function(seurat_obj, 
                             genes, 
                             filename, 
                             width = NULL, 
                             height = NULL, 
                             ncol = 3, 
                             mapping = NULL,
                             group_colors = NULL) {
  
  # 过滤存在于 Seurat 对象中的基因
  valid_genes <- intersect(genes, rownames(seurat_obj))
  
  if (length(valid_genes) == 0) {
    message(paste("No valid genes found for", filename))
    return(NULL)
  }
  
  # 计算行数和默认尺寸
  nrow <- ceiling(length(valid_genes) / ncol)
  
  # 创建绘图对象
  if (!is.null(mapping)) {
    p <- sc_violin(seurat_obj, 
                   features = valid_genes, 
                   .fun = function(d) dplyr::filter(d, value > 0), 
                   mapping = mapping, 
                   ncol = ncol) +
      guides(x = guide_axis(angle = -45)) +
      stat_compare_means(label = "p.signif")
    
    # 自动计算宽度
    if (is.null(width)) {
      x_var <- rlang::as_name(mapping$x)
      n_levels <- length(unique(seurat_obj@meta.data[[x_var]]))
      width <- 0.8 * n_levels * ncol
    }
  } else {
    p <- sc_violin(seurat_obj, 
                   features = valid_genes, 
                   .fun = function(d) dplyr::filter(d, value > 0),
                   ncol = ncol) +
      guides(x = guide_axis(angle = -45)) +
      stat_compare_means(label = "p.signif")
    
    # 默认宽度
    if (is.null(width)) {
      n_levels <- length(levels(Idents(seurat_obj)))
      width <- 0.8 * n_levels * ncol
    }
  }
  
  # 如果提供了group_colors，添加颜色标度
  if (!is.null(group_colors)) {
    p <- p + scale_fill_manual(values = group_colors)
  }
  
  # 自动计算高度
  if (is.null(height)) {
    height <- 4 * nrow
  }
  
  # 保存图片
  ggsave(paste0(filename, ".png"), plot = p, width = width, height = height, limitsize = FALSE)
  ggsave(paste0(filename, ".pdf"), plot = p, width = width, height = height, limitsize = FALSE)
  
  return(p)
}
# 使用示例：
## 节律基因 GO:0048511，基因太多了（300+）
# genes_CR <- c("Arntl", "Clock", "Cry1", "Cry2", "Csnk1d", "Csnk1e", "Npas2", "Nr1d1", "Per1", "Per2", "Per3")
## 铁硫簇组装基因 GO:0016226
# genes_FS <- c("Abcb7", "AK157302", "Bola2", "Bola3", "Ciao1", "Ciao2a", "Ciao2b", "Ciao3", 
#               "Ciapin1", "Fdx2", "Fxn", "Glrx3", "Glrx5", "Hscb", "Hspa9", "Iba57", "Isca1", 
#               "Isca2", "Iscu", "Lyrm4", "Ndor1", "Ndufab1", "Ndufab1-ps", "Nfs1", "Nfu1", 
#               "Nubp1", "Nubp2", "Nubpl", "Xdh")
## 定义颜色
# group_colors <- c(D = "#dd5633", L = "#5fb9d4")
## 定义 mapping
# mapping <- aes(x = ftype, y = value, fill = mfield)
#
## 单次绘图
# vlnplot_gene_set(seurat_object, genes = genes_CR, filename = "Vlnplot_CR", ncol = 4, 
#                  mapping = mapping, group_colors = group_colors)
#
## 批量绘图
## 定义基因集列表
# gene_sets <- list(
#   CR = list(genes = genes_CR, filename = "Vlnplot_CR", ncol = 4, mapping = mapping, group_colors = group_colors),
#   FS = list(genes = genes_FS, filename = "Vlnplot_FS", ncol = 6, mapping = mapping, group_colors = group_colors)
# )
## 批量绘图
# lapply(gene_sets, function(x) {
#   vlnplot_gene_set(seurat_object, 
#                    genes = x$genes, 
#                    filename = x$filename, 
#                    width = x$width, 
#                    height = x$height, 
#                    ncol = x$ncol, 
#                    mapping = x$mapping, 
#                    group_colors = group_colors)
# })


# 自定义颜色
custom_colors <- c("#5560AC", "#FCBB44", "#08306B", "#D3F0F2", "#67000d", 
                   "#ED6F6E", "#1F4527", "#D2D6F5", "#0D8B43", "#F4EEAC", 
                   "#FED98E", "#4758A2", "#E4E45F", "#C9DCC4", "#37939A", 
                   "#F28147", "#619CD9", "#EDADC5", "#F1766D", "#6CBEC3", 
                   "#ADAFB1", "#9ECAE1", "#68BD48", "#D8D9DA", "#7A70B5", 
                   "#ECB884", "#92A5D1", "#E08D8B", "#9584C1", "#AAD7C8", 
                   "#F1AEA7", "#7C9895", "#C5DFF4", "#304E7E", "#CEDBD2", 
                   "#AF8CBB", "#FAC074", "#9D9ECD", "#9FD4AE", "#D9B9D4")
