
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
plot_boxpair <- function(seurat_object, markers, target_genes, group_col = "group", save_plot = TRUE, 
                         filename = "boxpair", title = "", width = 4, height = 6) {
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
  expr_data <- FetchData(seurat_object, vars = c(valid_genes, group_col))
  
  # 4. 计算每个组内各目标基因的平均表达值，并转换为长格式数据
  avg_exp <- expr_data %>%
    group_by(across(all_of(group_col))) %>%
    summarise(across(all_of(valid_genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(cols = all_of(valid_genes),
                 names_to = "gene",
                 values_to = "expression")
  
  # 5. 绘制成对图
  p <- ggpaired(avg_exp, x = group_col, y = "expression", color = group_col, palette = "npg", point.size = 3) +
    stat_compare_means(aes_string(x = group_col, y = "expression"), paired = TRUE, method = "t.test", size = 4) +
    geom_text_repel(aes(label = gene), size = 4, box.padding = 0.5, point.padding = 0.3, segment.color = "grey50",
                    color = ifelse(avg_exp$gene %in% target_degs, "#8e44ad", "black")) +
    labs(title = title, x = "Group", y = "Average Expression", color = "Group") +
    theme_classic()
  
  # 显示图形
  print(p)
  
  # 6. 如果需要保存图形，则输出为 PDF 和 PNG 文件
  if (save_plot) {
    ggsave(filename = paste0(filename, ".pdf"), plot = p, width = width, height = height)
    ggsave(filename = paste0(filename, ".png"), plot = p, width = width, height = height)
  }
  
  # 返回绘图对象，方便后续调整
  return(p)
}
# 用法示例：
# plot_boxpair(seurat_object, markers = markers, target_genes = genes_CR,
#              group_col = "group", save_plot = T, filename = "boxpair_CR", width = 6, height = 8)


# 绘制 UMAP/T-SNE 的函数
plot_seurat_dim <- function(seurat_obj, 
                            reduction = "umap", 
                            ident = "orig.ident", 
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
  p <- sc_dim(seurat_obj, reduction = reduction) + 
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
# plot_seurat_dim(seurat_object, ident = "RNA_snn_res.0.8", filename = "UMAP_Cluster")
# 优先级顺序：name_to_color > custom_colors >  viridis 调色板
# name_to_color = c("Fibroblasts" = "blue", "Other Cells" = "#FF5733")
# plot_seurat_dim(seurat_object, reduction = "umap", ident = "RNA_snn_res.0.8", 
#                 height = 6, width = 7, filename = "UMAP_Cluster",
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


# 自定义颜色
custom_colors <- c("#5560AC", "#FCBB44", "#08306B", "#D3F0F2", "#67000d", 
                   "#ED6F6E", "#1F4527", "#D2D6F5", "#0D8B43", "#F4EEAC", 
                   "#304E7E", "#FED98E", "#4758A2", "#E4E45F", "#C9DCC4", 
                   "#37939A", "#F28147", "#619CD9", "#EDADC5", "#F1766D", 
                   "#6CBEC3", "#ADAFB1", "#9ECAE1", "#68BD48", "#D8D9DA", 
                   "#7A70B5", "#ECB884", "#92A5D1", "#E08D8B", "#9584C1", 
                   "#AAD7C8", "#F1AEA7", "#7C9895", "#C5DFF4", "#CEDBD2", 
                   "#AF8CBB", "#FAC074", "#9D9ECD", "#9FD4AE", "#D9B9D4")
