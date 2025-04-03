
############################################## 单细胞分析函数 ############################################

#' 执行 Seurat 对象的初步质量控制
#'
#' 该函数计算线粒体基因、核糖体基因和血红蛋白基因的百分比，
#' 并生成小提琴图展示这些QC指标。
#'
#' @param seurat_obj 包含单细胞RNA-seq数据的Seurat对象
#' @param dir 输出图片的保存目录路径。如果为NULL，则不保存图片
#'
#' @return 返回添加了QC指标的Seurat对象，包含以下新列：
#' \itemize{
#'   \item percent_mito - 线粒体基因百分比
#'   \item percent_ribo - 核糖体基因百分比 
#'   \item percent_hb - 血红蛋白基因百分比
#' }
#'
#' @details
#' 该函数会：
#' \enumerate{
#'   \item 计算线粒体基因(以"mt-"开头)的百分比
#'   \item 计算核糖体基因(以"Rp[sl]"开头)的百分比
#'   \item 计算血红蛋白基因(以"Hb"开头但不包含"p")的百分比
#'   \item 生成包含nFeature_RNA、nCount_RNA和上述百分比的小提琴图
#' }
#'
#' @examples
#' \dontrun{
#' seurat_object <- basic_qc(seurat_object)
#' seurat_object <- basic_qc(seurat_object, dir = "QC_results")
#' }
#'
#' @export
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


#' 绘制组间基因表达差异的配对箱线图
#'
#' 该函数用于可视化不同组别间基因表达差异，特别适用于配对样本分析，
#' 会自动标注差异显著的基因并执行配对t检验。
#'
#' @param seurat_object 包含单细胞 RNA-seq 数据的 Seurat 对象
#' @param markers 差异表达分析结果数据框，需包含 p_val_adj 列
#' @param target_genes 需要可视化的目标基因列表
#' @param group 指定分组变量的列名（默认："group"）
#' @param save_plot 是否保存图片（默认：TRUE）
#' @param filename 输出文件名前缀（默认："boxpair"）
#' @param title 图片标题（默认：空）
#' @param width 图片宽度（单位：英寸，默认：4）
#' @param height 图片高度（单位：英寸，默认：6）
#' @param group_colors 自定义分组颜色向量（默认：NULL，使用 npg 调色板）
#'
#' @return 返回一个 ggplot 对象，包含以下元素：
#' \itemize{
#'   \item 配对箱线图展示各组基因平均表达
#'   \item 显著差异基因用紫色标注（p_val_adj < 0.05）
#'   \item 自动添加配对t检验p值
#' }
#'
#' @details
#' 函数执行流程：
#' \enumerate{
#'   \item 筛选存在于 seurat_object 中的目标基因
#'   \item 筛选在 markers 中显著差异的基因（p_val_adj < 0.05）
#'   \item 提取表达矩阵和分组信息
#'   \item 计算各组基因平均表达值
#'   \item 绘制配对箱线图并添加统计检验
#'   \item 可选保存为 PDF 和 PNG 格式
#' }
#'
#' @examples
#' \dontrun{
#' # 基本用法
#' plot_boxpair(seurat_obj, markers = deg_results, 
#'              target_genes = c("CD4", "CD8A", "FOXP3"))
#'              
#' # 自定义参数
#' plot_boxpair(seurat_obj, markers = deg_results, 
#'              target_genes = top_genes, group = "treatment",
#'              group_colors = c("Control" = "grey", "Treatment" = "red"))
#' }
#'
#' @importFrom dplyr filter group_by summarise across
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggpaired stat_compare_means
#' @importFrom ggrepel geom_text_repel
#' @export
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


#' 绘制基因集的表达特征图
#'
#' 该函数用于可视化一组基因在单细胞数据中的表达模式，自动调整图片尺寸并保存为多种格式。
#'
#' @param seurat_obj 包含 scRNA-seq 数据的 Seurat 对象
#' @param genes 需要可视化的基因列表
#' @param output_dir 输出目录路径
#' @param file_prefix 输出文件名前缀（用于区分不同基因集）
#' @param base_width 单个基因子图的基准宽度（单位：英寸，默认：2）
#' @param base_height 单个基因子图的基准高度（单位：英寸，默认：2）
#' @param ncol 每行显示的基因数量（默认：3）
#'
#' @return 无返回值，直接输出图片文件到指定目录
#'
#' @details
#' 函数执行流程：
#' \enumerate{
#'   \item 检查并创建输出目录（如不存在）
#'   \item 筛选存在于 seurat_obj 中的有效基因
#'   \item 根据基因数量和布局参数自动计算图片尺寸
#'   \item 生成特征表达图并添加标题
#'   \item 保存为 PNG 和 PDF 格式（文件名自动添加前缀）
#' }
#'
#' @note
#' 重要特性：
#' \itemize{
#'   \item 自动跳过不存在的基因并提示
#'   \item 图片尺寸根据基因数量和布局参数动态调整
#'   \item 关闭了 ggsave 的尺寸限制（limitsize = FALSE）
#' }
#'
#' @examples
#' \dontrun{
#' # 基本用法
#' featureplot_gene_set(seurat_obj, 
#'                     genes = c("CD4", "CD8A", "FOXP3"),
#'                     output_dir = "results/plots",
#'                     file_prefix = "Tcell_markers")
#'                     
#' # 自定义布局
#' featureplot_gene_set(seurat_obj, 
#'                     genes = c("CD4", "CD8A", "FOXP3", "IL2RA"),
#'                     output_dir = "results/plots",
#'                     file_prefix = "Tcell_markers",
#'                     ncol = 2, 
#'                     base_width = 3)
#' }
#'
#' @importFrom Seurat FeaturePlot
#' @export
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


#' 绘制 Seurat 对象的降维可视化图（UMAP/t-SNE）
#'
#' 该函数提供灵活的降维可视化功能，支持自定义颜色、标签和细胞数量显示。
#'
#' @param seurat_obj 包含单细胞 RNA-seq 数据的 Seurat 对象
#' @param reduction 使用的降维方法（"umap"或"tsne"，默认："umap"）
#' @param ident 用于分组的元数据列名（默认："orig.ident"）
#' @param pointsize 点的大小（默认：1）
#' @param height 图片高度（单位：英寸，默认：6）
#' @param width 图片宽度（单位：英寸，默认：7）
#' @param filename 输出文件名前缀（默认："Dimplot"）
#' @param custom_colors 自定义颜色向量（长度需≥分组数，默认：NULL）
#' @param name_to_color 命名颜色向量（直接为每个分组指定颜色，默认：NULL）
#' @param label_color 标签文本颜色（默认：'#5D478B'）
#' @param legend_size 图例点的大小（默认：5）
#' @param color_palette viridis 调色板选项（默认：'H'）
#' @param show_cellnum 是否显示各分组的细胞数量（默认：FALSE）
#' @param title 图片标题（默认：""）
#' @param legend_title 图例标题（默认使用 ident 参数值）
#' @param save_plot 是否保存图片（默认：TRUE）
#'
#' @return 返回包含以下元素的列表：
#' \itemize{
#'   \item plot - ggplot 对象
#'   \item color_map - 使用的颜色映射关系
#' }
#'
#' @details
#' 颜色设置优先级：
#' \enumerate{
#'   \item name_to_color（直接指定分组颜色）
#'   \item custom_colors（自定义颜色向量）
#'   \item viridis 调色板（默认）
#' }
#'
#' @section 功能特性：
#' \itemize{
#'   \item 自动调整图片宽度当显示细胞数量时
#'   \item 支持三种颜色设置方式
#'   \item 可选的细胞数量标注
#'   \item 返回颜色映射关系便于后续使用
#' }
#'
#' @examples
#' \dontrun{
#' # 基本用法
#' result <- plot_seurat_dim(seurat_object, ident = "RNA_snn_res.0.8")
#' 
#' # 显示细胞数量
#' plot_seurat_dim(seurat_object, show_cellnum = TRUE)
#' 
#' # 自定义颜色
#' plot_seurat_dim(seurat_object,
#'                name_to_color = c("Fibroblasts" = "blue", "Other Cells" = "#FF5733"))
#' 
#' # 使用 tsne 降维
#' plot_seurat_dim(seurat_object, reduction = "tsne")
#' 
#' plot_seurat_dim(seurat_object, ident = "RNA_snn_res.0.8", filename = "UMAP_cluster")
#' 
#' name_to_color = c("Fibroblasts" = "blue", "Other Cells" = "#FF5733")
#' plot_seurat_dim(seurat_object, reduction = "umap", ident = "RNA_snn_res.0.8", 
#'                 height = 6, width = 7, filename = "UMAP_cluster",
#'                 custom_colors = NULL, name_to_color = NULL, 
#'                 label_color = '#5D478B', legend_size = 5, color_palette = 'H',
#'                 title = "", legend_title = ident)
#' }
#'
#' @importFrom ggsc sc_dim sc_dim_geom_label
#' @importFrom ggrepel geom_text_repel
#' @importFrom viridis viridis_pal
#' @importFrom ggplot2 scale_color_manual guide_legend guides
#' @export
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
                            show_cellnum = FALSE,
                            title = "",
                            legend_title = ident, 
                            save_plot = TRUE) { 
  library(ggsc)
  library(ggrepel)
  
  # 设置身份标识符
  Idents(seurat_obj) <- ident
  ids <- Idents(seurat_obj)
  unique_ids <- levels(ids)
  
  # 如果 show_cellnum = TRUE，计算各分组的细胞数量
  if (show_cellnum) {
    cell_counts <- table(ids)[unique_ids]  # 按排序后的顺序提取
    new_labels <- paste0(unique_ids, " (", cell_counts, ")")
    width <- width + 1
  } else {
    new_labels <- unique_ids
  }
  
  # 自动颜色映射逻辑
  if (!is.null(name_to_color)) {
    color_scale <- scale_color_manual(values = name_to_color, labels = new_labels)
    color_map <- as.list(name_to_color)
  } else if (!is.null(custom_colors)) {
    if (length(unique_ids) > length(custom_colors)) {
      stop("颜色数量少于簇的数量！")
    }
    color_map <- setNames(custom_colors, unique_ids)
    color_scale <- scale_color_manual(values = color_map, labels = new_labels)
  } else {
    color_map <- viridis::viridis_pal(option = color_palette)(length(unique_ids))
    names(color_map) <- unique_ids
    color_scale <- scale_color_manual(values = color_map, labels = new_labels)
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


#' 绘制细胞类型标记基因表达点图
#'
#' 该函数用于可视化不同细胞类型的标记基因在各分组中的表达模式，自动调整图片尺寸并保存为 PDF 文件。
#'
#' @param seurat_obj 包含 scRNA-seq 数据的 Seurat 对象
#' @param marker_list 包含细胞类型标记基因的列表，格式见示例
#' @param group.by 用于分组的元数据列名
#' @param output_dir 输出目录路径（默认："."）
#'
#' @return 返回包含所有绘图对象的列表（不可见），同时自动保存 PDF 文件到指定目录
#'
#' @details
#' 函数执行流程：
#' \enumerate{
#'   \item 检查每个细胞类型的标记基因是否存在
#'   \item 对每个细胞类型生成点图（Dotplot）
#'   \item 自动调整图片尺寸（基于基因数量）
#'   \item 保存为 PDF 文件（文件名包含细胞类型信息）
#' }
#'
#' @section 图片特性：
#' \itemize{
#'   \item 点的大小表示表达比例
#'   \item 颜色深浅表示平均表达量（绿色到红色渐变）
#'   \item X 轴标签 45 度倾斜提高可读性
#'   \item 标题加粗显示细胞类型
#' }
#'
#' @section 动态尺寸调整规则：
#' \itemize{
#'   \item 宽度：max(3, n_genes * 1.2) 英寸
#'   \item 高度：max(4, n_genes * 0.3) 英寸
#'   \item 关闭 ggsave 的尺寸限制（limitsize = FALSE）
#' }
#'
#' @examples
#' \dontrun{
#' # 准备标记基因列表
#' marker_genes <- list(
#'   T_Cell = c("Cd3g","Cd3e","Cd4","Cd8a"),
#'   B_Cell = c("Cd19", "Cd79a", "Cd79b"),
#'   Macrophage = c("Adgre1","Csf1r","Trem2")
#' )
#'
#' # 基本用法
#' plot_markers(seurat_object, 
#'             marker_list = marker_genes,
#'             group.by = "RNA_snn_res.0.8")
#'
#' # 指定输出目录
#' plot_markers(seurat_object,
#'             marker_genes,
#'             group.by = "celltype",
#'             output_dir = "results/marker_plots")
#' }
#'
#' @importFrom ggsc sc_dot
#' @importFrom ggplot2 scale_color_gradient ggtitle theme element_text
#' @export
plot_markers <- function(seurat_obj, marker_list, group.by, output_dir = ".") {
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
      filename = paste0(output_dir, "/check_", gsub("[/]", "_", cell_type), ".pdf"),
      plot = p,
      width = base_width,
      height = base_height,
      limitsize = FALSE
    )
    
    return(p)
  })
}


#' 绘制基因集表达小提琴图
#'
#' 该函数用于可视化一组基因在不同分组中的表达分布，支持自定义分组映射和颜色，自动调整图片尺寸。
#'
#' @param seurat_obj 包含 scRNA-seq 数据的 Seurat 对象
#' @param genes 需要可视化的基因向量
#' @param filename 输出文件名前缀（不含扩展名）
#' @param width 图片宽度（单位：英寸，NULL 时自动计算）
#' @param height 图片高度（单位：英寸，NULL 时自动计算）
#' @param ncol 每行显示的子图数量（默认：3）
#' @param mapping ggplot2 映射关系（aes 对象），用于自定义 x 轴和填充变量
#' @param group_colors 自定义分组颜色向量
#'
#' @return 返回 ggplot 对象，同时自动保存 PNG 和 PDF 文件
#'
#' @details
#' 核心功能：
#' \itemize{
#'   \item 自动过滤不存在的基因
#'   \item 支持两种绘图模式：默认分组或自定义映射
#'   \item 自动添加组间显著性比较（p.signif）
#'   \item X 轴标签 45 度倾斜提高可读性
#' }
#'
#' @section 自动尺寸计算规则：
#' \itemize{
#'   \item 宽度：0.8 * 分组数量 * ncol（当 width = NULL 时）
#'   \item 高度：4 * 行数（当 height = NULL 时）
#'   \item 行数：ceiling(基因数量 / ncol)
#' }
#'
#' @section 高级用法：
#' \itemize{
#'   \item 使用 mapping 参数可自定义分组变量和填充变量
#'   \item 通过 group_colors 自定义颜色方案
#'   \item 支持批量处理多个基因集（见示例）
#' }
#'
#' @examples
#' \dontrun{
#' # 基本用法（使用默认分组）
#' vlnplot_gene_set(seurat_obj, 
#'                 genes = c("Gene1", "Gene2"), 
#'                 filename = "my_genes")
#'
#' # 自定义分组和颜色
#' my_mapping <- aes(x = cell_type, y = value, fill = condition)
#' my_colors <- c("Type1" = "#1f77b4", "Type2" = "#ff7f0e")
#' 
#' vlnplot_gene_set(seurat_obj,
#'                 genes = marker_genes,
#'                 filename = "custom_vln",
#'                 mapping = my_mapping,
#'                 group_colors = my_colors)
#'
#' # 批量处理多个基因集
#' gene_sets <- list(
#'   set1 = list(genes = c("Gene1", "Gene2"), filename = "set1"),
#'   set2 = list(genes = c("Gene3", "Gene4"), filename = "set2")
#' )
#' 
#' lapply(gene_sets, function(x) {
#'   vlnplot_gene_set(seurat_obj, 
#'                   genes = x$genes, 
#'                   filename = x$filename)
#' })
#' }
#'
#' @importFrom ggsc sc_violin
#' @importFrom ggplot2 aes guide_axis scale_fill_manual
#' @importFrom ggpubr stat_compare_means
#' @importFrom dplyr filter
#' @export
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


# 自定义颜色
custom_colors <- c("#5560AC", "#FCBB44", "#08306B", "#D3F0F2", "#67000d", 
                   "#ED6F6E", "#1F4527", "#D2D6F5", "#0D8B43", "#F4EEAC", 
                   "#FED98E", "#4758A2", "#E4E45F", "#C9DCC4", "#37939A", 
                   "#F28147", "#619CD9", "#EDADC5", "#F1766D", "#6CBEC3", 
                   "#ADAFB1", "#9ECAE1", "#68BD48", "#D8D9DA", "#7A70B5", 
                   "#ECB884", "#92A5D1", "#E08D8B", "#9584C1", "#AAD7C8", 
                   "#F1AEA7", "#7C9895", "#C5DFF4", "#304E7E", "#CEDBD2", 
                   "#AF8CBB", "#FAC074", "#9D9ECD", "#9FD4AE", "#D9B9D4")
