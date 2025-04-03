
library(Seurat)
library(COSG)
library(ggsc)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggalluvial)
library(openxlsx)

# 清理环境
rm(list = ls())
gc()

# ======================================= 参数设置 =======================================
# 基础目录设置
base_dir <- "~/scLC"
project_name <- "W_Penis"  # 项目名称

seurat_obj <- "seurat_W_Penis.rds"
clustering_resolution <- "0.5"  # 聚类分辨率

# 子目录自动生成
datadir <- file.path(base_dir, "data")
resultdir <- file.path(base_dir, "results")
scriptdir <- file.path(base_dir, "scripts")
project_dir <- file.path(resultdir, project_name)
analysis_dir <- file.path(project_dir, paste0("check-markers-", clustering_resolution))
checkmarker_dir <- file.path(analysis_dir, "checkplot")

# 用于细胞注释的文件（需手动整理，参考示例见文件夹）
marker_file <- file.path(datadir, "markers_penis.csv")  # 细胞标记基因文件，必须含有 celltype、ftype、markers 列
mapping_file <- file.path(analysis_dir, "cluster_name.csv")  # 细胞注释映射文件，必须含有 cluster、celltype、ftype 列

# 创建所需目录
for (dir in c(analysis_dir, checkmarker_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# FindAllMarkers 参数
analysis_params <- list(
  min_pct = 0.25,
  logfc_threshold = 0.25,
  top_n_genes = 100,
  p_val_adj_cutoff = 0.05
)

# 绘图参数
plot_params <- list(
  point_size = 1,
  umap_width = 8,
  umap_height = 6,
  cell_palette = NULL,  # 设置为 NULL 时默认使用色盲友好颜色
  group_colors = c("D" = "#dd5633", "L" = "#5fb9d4"),
  cluster_var = paste0("RNA_snn_res.", clustering_resolution),
  celltype_var = "celltype",
  parent_type_var = "ftype",
  condition_var = "mfield"
)


# ================================= 以下部分都不需要修改 =================================
# ========================== 但需要在注释前整理细胞注释映射文件 ==========================
# ======================================= 函数引入 =======================================
source(file.path(scriptdir, "Functions_scRNA.R"))

# ======================================= 辅助函数 =======================================
# 计算差异基因并保存为 Excel 的函数
calculate_and_save_markers <- function(seurat_obj, ident_var, params) {
  Idents(seurat_obj) <- ident_var
  
  # 计算差异基因
  message("计算差异基因...")
  markers <- FindAllMarkers(
    seurat_obj, 
    only.pos = TRUE, 
    min.pct = params$min_pct, 
    logfc.threshold = params$logfc_threshold
  )
  
  # 筛选并排序
  markers_top <- markers %>%
    filter(p_val_adj < params$p_val_adj_cutoff) %>%
    arrange(cluster, desc(abs(avg_log2FC))) %>% 
    group_by(cluster) %>% 
    slice_head(n = params$top_n_genes) %>% 
    ungroup()
  
  # 创建 Excel 工作簿
  wb <- createWorkbook()
  clusters <- unique(markers_top$cluster)
  
  # 按 cluster 分 sheet 保存
  for (clust in clusters) {
    sheet_data <- markers_top %>% filter(cluster == clust)
    addWorksheet(wb, sheetName = as.character(clust))
    writeData(wb, sheet = as.character(clust), x = sheet_data)
  }
  
  # 保存工作簿
  output_file <- paste0("DEGs_clusters_top", params$top_n_genes, ".xlsx")
  saveWorkbook(wb, file = output_file, overwrite = TRUE)
  message(paste("差异基因已保存至:", output_file))
  
  return(markers_top)
}

# 提取并保存 UMAP 坐标
extract_and_save_umap <- function(seurat_obj, cluster_var) {
  # 提取坐标
  umap <- seurat_obj@reductions$umap@cell.embeddings %>%
    as.data.frame() %>% 
    cbind(cluster = seurat_obj@meta.data[[cluster_var]]) %>% 
    rownames_to_column("cellid")
  
  # 保存全部坐标
  write.csv(umap, "umap_all.csv", row.names = FALSE)
  
  # 计算并保存每个 cluster 的平均坐标
  umap_means <- aggregate(cbind(umap_1, umap_2) ~ cluster, data = umap, FUN = mean)
  write.csv(umap_means, "umap_cluster.csv", row.names = FALSE)
  
  return(list(umap = umap, umap_means = umap_means))
}

# 细胞类型标记基因处理函数
process_cell_markers <- function(marker_file) {
  # 读取细胞标记基因
  anno_markers_df <- read.csv(marker_file)
  
  # 处理标记基因列表
  anno_markers <- setNames(
    lapply(anno_markers_df$markers, function(x) unlist(strsplit(x, ", "))),
    gsub("[/ ]", "_", anno_markers_df$subtype)
  )
  
  return(list(anno_markers_df = anno_markers_df, anno_markers = anno_markers))
}

# 计算 markers 与 cluster 差异基因的重叠度
calculate_marker_overlap <- function(marker_gene_list, cluster_gene_list) {
  # 计算交集矩阵
  overlap_matrix <- sapply(cluster_gene_list, function(cluster_genes) {
    sapply(marker_gene_list, function(celltype_genes) {
      length(intersect(tolower(cluster_genes), tolower(celltype_genes)))
    })
  })
  
  # 保存为数据框
  overlap_df <- as.data.frame(overlap_matrix)
  
  # 绘制热图
  pdf("Cluster_CellType_Marker_Overlap.pdf", width = 9, height = 7)
  pheatmap(overlap_matrix,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           display_numbers = TRUE,
           main = "Cluster-CellType Marker Overlap")
  dev.off()
  
  return(overlap_df)
}

# 绘制所有标记基因的点图
plot_all_markers <- function(seurat_obj, anno_markers, group_by_var) {
  all_check_genes <- unique(unlist(anno_markers))
  valid_genes <- all_check_genes[all_check_genes %in% rownames(seurat_obj)]
  
  p <- sc_dot(
    object = seurat_obj,
    features = valid_genes,
    group.by = group_by_var, 
    dot.scale = 6
  ) + coord_flip() + 
    scale_color_gradient(low = "#4d8076", high = "red")
  
  # 动态调整图形大小
  h <- length(valid_genes) / 6 + 3
  w <- length(unique(seurat_obj@meta.data[[group_by_var]])) * 0.3 + 1
  
  ggsave(filename = "check_for_all.pdf", plot = p, width = w, height = h)
  
  return(p)
}

# 基于细胞类型注释映射表注释 cells
annotate_cells <- function(seurat_obj, mapping_file, cluster_var) {
  # 从文件加载映射
  cluster_name <- read.csv(mapping_file)
  cmeta <- seurat_obj@meta.data
  
  # 映射细胞类型
  cmeta$celltype <- cluster_name$celltype[match(cmeta[[cluster_var]], cluster_name$cluster)] %>% as.factor()
  cmeta$ftype <- cluster_name$ftype[match(cmeta[[cluster_var]], cluster_name$cluster)] %>% as.factor()
  
  # 更新 metadata
  seurat_obj@meta.data <- cmeta

  return(seurat_obj)
}

# 绘制不同注释的 UMAP 图
plot_umap_for_group <- function(plotgroup) {
  num <- length(levels(seurat_object@meta.data[[plotgroup]]))
  result <- plot_seurat_dim(seurat_object, 
                            ident = plotgroup, 
                            pointsize = 1, 
                            custom_colors = custom_colors[1:num], 
                            filename = paste0("UMAP_", plotgroup), 
                            legend_title = plotgroup,
                            width = 8, height = 6)
  
  return(result$color_map)
}

# 比较不同组间的细胞类型分布
compare_groups <- function(seurat_obj, group_var, condition_var, color_map) {
  # 按条件分组的 UMAP
  mapping <- aes(color = !!sym(group_var))
  p <- sc_dim(seurat_obj, reduction = "umap", mapping = mapping, pointsize = 1.2) +
    facet_wrap(~ get(condition_var), scales = "free") +
    sc_dim_geom_label(geom = ggrepel::geom_text_repel, color = '#5D478B', bg.color = 'white') + 
    scale_color_manual(values = color_map) + 
    guides(color = guide_legend(
      title = "Cell Type", 
      override.aes = list(size = 5))) +
    ggtitle(paste("UMAP by", condition_var))
  
  ggsave(paste0("Diff_group_", group_var, ".pdf"), plot = p, width = 13, height = 6)
  ggsave(paste0("Diff_group_", group_var, ".png"), plot = p, width = 13, height = 6)
  
  # 计算各组不同细胞群比例
  cellratio <- prop.table(table(seurat_obj@meta.data[[group_var]], 
                                seurat_obj@meta.data[[condition_var]]), margin = 2) %>% 
    as.data.frame()
  names(cellratio) <- c(group_var, condition_var, "Freq")
  
  # 绘制堆叠柱状图
  p <- ggplot(cellratio, aes(x = .data[[condition_var]], y = Freq, 
                             fill = .data[[group_var]], 
                             stratum = .data[[group_var]], 
                             alluvium = .data[[group_var]])) + 
    geom_col(width = 0.6, color = NA) + 
    geom_flow(width = 0.5, alpha = 0.5) + 
    theme_minimal() + 
    labs(x = condition_var, y = "Proportion", title = "Cell Type Ratio Across Groups") + 
    scale_fill_manual(values = color_map)
  
  ggsave(paste0("percentage_", group_var, ".pdf"), plot = p, width = 6, height = 6)
  ggsave(paste0("percentage_", group_var, ".png"), plot = p, width = 6, height = 6)
  
  return(list(plot = p, data = cellratio))
}

# ======================================= 数据加载 =======================================
setwd(analysis_dir)
seurat_object <- readRDS(file.path(project_dir, seurat_obj))

# ==================================== 手动注释前准备 ====================================
# 计算并保存差异基因
markers_top <- calculate_and_save_markers(
  seurat_obj = seurat_object, 
  ident_var = plot_params$cluster_var, 
  params = analysis_params
)

# 提取并保存 UMAP 坐标
extract_and_save_umap(seurat_obj = seurat_object, cluster_var = plot_params$cluster_var)

# 处理细胞类型标记基因
cell_markers <- process_cell_markers(marker_file)

# 计算标记基因与 cluster 的重叠度
cluster_genes <- split(markers_top$gene, markers_top$cluster)
overlap_df <- calculate_marker_overlap(cell_markers$anno_markers, cluster_genes)

# 绘制标记基因点图
message("绘制标记基因点图...")
plot_markers(seurat_object, cell_markers$anno_markers, group.by = plot_params$cluster_var, output_dir = checkmarker_dir)
plot_all_markers(seurat_obj = seurat_object, anno_markers = cell_markers$anno_markers, group_by_var = plot_params$cluster_var)

# 使用 cluster 注释映射文件
message("使用映射文件进行细胞类型注释...")
seurat_object <- annotate_cells(seurat_obj = seurat_object, 
                                mapping_file = mapping_file, 
                                cluster_var = plot_params$cluster_var)


# 要绘制的分组
plotgroups <- c(plot_params$cluster_var, plot_params$celltype_var, plot_params$parent_type_var)
color_maps <- lapply(plotgroups, plot_umap_for_group)
names(color_maps) <- plotgroups  # 命名列表

# 比较不同组之间的细胞分布
message("比较不同组间的细胞类型分布...")
Idents(seurat_object) <- seurat_object@meta.data[[plot_params$condition_var]]
comparison_results <- lapply(plotgroups, function(group_var) {
  compare_groups(
    seurat_obj = seurat_object,  
    group_var = group_var,
    condition_var = plot_params$condition_var,
    color_map = color_maps[[group_var]]
  )
})
