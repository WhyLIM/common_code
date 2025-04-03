
library(Seurat)
library(ggsc)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggtangle)
library(clusterProfiler)
library(rlang)

# 清理环境
rm(list = ls())
gc()

# ======================================= 参数设置 =======================================
# 基础目录设置
base_dir <- "~/scLC"
project_name <- "W_Penis"  # 项目名称
clustering_resolution <- "0.5"  # 聚类分辨率
seurat_obj <- "seurat_celltype.rds"  # Seurat 对象

# 子目录自动生成
datadir <- file.path(base_dir, "data")
resultdir <- file.path(base_dir, "results")
scriptdir <- file.path(base_dir, "scripts")
project_dir <- file.path(resultdir, project_name)
analysis_dir <- file.path(project_dir, paste0("check-markers-", clustering_resolution))
keygenes_dir <- file.path(analysis_dir, "KeyGenes")

# ====================================== 基因集定义 ======================================
# 只需定义基因集名称和基因列表，其他参数会自动生成
gene_sets_raw <- list(
  # 节律基因 (Circadian Rhythm)
  CR = c("Arntl", "Clock", "Cry1", "Cry2", "Csnk1d", "Csnk1e", "Npas2", "Nr1d1", "Per1", "Per2", "Per3"),
  
  # 铁硫簇组装基因 (Fe-S Cluster)
  FS = c("Abcb7", "AK157302", "Bola2", "Bola3", "Ciao1", "Ciao2a", "Ciao2b", "Ciao3", 
         "Ciapin1", "Fdx2", "Fxn", "Glrx3", "Glrx5", "Hscb", "Hspa9", "Iba57", "Isca1", 
         "Isca2", "Iscu", "Lyrm4", "Ndor1", "Ndufab1", "Ndufab1-ps", "Nfs1", "Nfu1", 
         "Nubp1", "Nubp2", "Nubpl", "Xdh"),
  
  # 前列腺素合成基因（Prostaglandin synthesis）
  SP = c("Anxa1", "Avp", "Avpr1a", "Cd74", "Cthrc1", "Daglb", "Edn1", "Edn2", "Fabp5", 
         "Hpgds", "Il1b", "Mapk9", "Mif", "Pibf1", "Pla2g2a", "Pla2g3", "Pla2g4a", 
         "Pla2g4f", "Pla2g10", "Pnpla8", "Prxl2b", "Ptgds", "Ptges", "Ptges2", 
         "Ptges3", "Ptges3-ps", "Ptgis", "Ptgs1", "Ptgs2", "Sco1", "Sirt1", "Sphk1", 
         "Tbxas1"),
  
  # 精囊分泌蛋白基因（Seminal vesicle secretion protein）
  SV = c("Svs4", "Svs3b", "Aoc1l3", "Semg1", "Svs6", "Pate4", "Svs3a", "Svs5")
  
  # 添加新基因集只需添加名称和基因列表
  # NEW = c("Gene1", "Gene2", "Gene3")
)

# 自动生成完整的基因集配置
create_gene_sets <- function(gene_sets_raw) {
  gene_sets <- list()
  
  for (set_name in names(gene_sets_raw)) {
    # 根据基因数量确定每行显示的列数
    gene_count <- length(gene_sets_raw[[set_name]])
    ncols <- if (gene_count <= 12) 4 else if (gene_count <= 24) 5 else 6
    
    # 创建完整的基因集配置
    gene_sets[[set_name]] <- list(
      genes = gene_sets_raw[[set_name]],
      vln_filename = paste0("Vlnplot_", set_name),
      feature_prefix = set_name,
      boxpair_dir = paste0("boxpair_plots/", set_name),
      ncol = ncols
    )
  }
  
  return(gene_sets)
}

# 生成基因集配置
gene_sets <- create_gene_sets(gene_sets_raw)
# 自动创建需要的目录
create_directories <- function(keygenes_dir, gene_sets) {
  # 创建主目录
  if (!dir.exists(keygenes_dir)) {
    dir.create(keygenes_dir, recursive = TRUE)
  }
  
  # 创建特征图目录
  featureplots_dir <- file.path(keygenes_dir, "featureplots")
  if (!dir.exists(featureplots_dir)) {
    dir.create(featureplots_dir, recursive = TRUE)
  }
  
  # 为每个基因集创建箱线图目录
  for (set_name in names(gene_sets)) {
    set_info <- gene_sets[[set_name]]
    boxpair_path <- file.path(keygenes_dir, set_info$boxpair_dir)
    if (!dir.exists(boxpair_path)) {
      dir.create(boxpair_path, recursive = TRUE)
    }
  }
}

# 创建所需目录
create_directories(keygenes_dir, gene_sets)

# ======================================= 绘图配置 =======================================
# 绘图参数
plot_params <- list(
  group_colors = c(D = "#dd5633", L = "#5fb9d4"),
  group_var = "ftype",
  condition_var = "mfield"
)


# ================================= 以下部分都不需要修改 =================================
# ======================================= 函数引入 =======================================
source(file.path(scriptdir, "Functions_scRNA.R"))
source(file.path(scriptdir, "Functions_Enrich.R"))

# ======================================= 加载数据 =======================================
# 读取 Seurat 对象
seurat_object <- readRDS(file.path(analysis_dir, seurat_obj))

setwd(keygenes_dir)

# ======================================= 小提琴图 =======================================
# 基础映射
base_mapping <- aes(x = !!sym(plot_params$group_var), y = value, fill = !!sym(plot_params$condition_var))

# 批量绘制小提琴图
for (set_name in names(gene_sets)) {
  set_info <- gene_sets[[set_name]]
  vlnplot_gene_set(
    seurat_object, 
    genes = set_info$genes, 
    filename = set_info$vln_filename, 
    ncol = set_info$ncol, 
    mapping = base_mapping, 
    group_colors = plot_params$group_colors
  )
}

# ====================================== 配对箱线图 ======================================
# 根据 group 变量分组
Idents(seurat_object) <- plot_params$group_var
celltypes <- levels(seurat_object@meta.data[[plot_params$group_var]])
seurat_cells <- list()

for (ct in celltypes) {
  seurat_cells[[ct]] <- subset(seurat_object, subset = !!sym(plot_params$group_var) == ct)
}

# 读取差异基因
deg_list <- readRDS(file.path(analysis_dir, "DEA/deg_list.rds"))

# 批量绘制配对箱线图
for (set_name in names(gene_sets)) {
  set_info <- gene_sets[[set_name]]
  boxpair_gene_set(
    seurat_list = seurat_cells,
    deg_list = deg_list,
    genes = set_info$genes,
    output_dir = set_info$boxpair_dir,
    group = plot_params$condition_var,
    group_colors = plot_params$group_colors,
    base_height = 0.5,
    base_width = 6
  )
}

# ====================================== 基因表达图 ======================================
# 批量绘制特征图
for (set_name in names(gene_sets)) {
  set_info <- gene_sets[[set_name]]
  featureplot_gene_set(
    seurat_obj = seurat_object,
    genes = set_info$genes,
    output_dir = "featureplots",
    file_prefix = set_info$feature_prefix,
    ncol = set_info$ncol
  )
}
