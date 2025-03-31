
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

# ======================================= 参数设置（可自定义） =======================================
# 基础目录设置 - 只需修改这个根目录，其他目录会自动生成
base_dir <- "~/scLC"
project_name <- "W_Heart"  # 项目名称，可以修改
clustering_resolution <- "0.8"  # 聚类分辨率，可以修改

# 子目录自动生成
datadir <- file.path(base_dir, "data")
resultdir <- file.path(base_dir, "results")
scriptdir <- file.path(base_dir, "scripts")
project_dir <- file.path(resultdir, project_name)
analysis_dir <- file.path(project_dir, paste0("check-markers-", clustering_resolution))
keygenes_dir <- file.path(analysis_dir, "KeyGenes")

# 自动创建需要的目录
for (dir in c(keygenes_dir, 
              file.path(keygenes_dir, "boxpair_plots", "CR"),
              file.path(keygenes_dir, "boxpair_plots", "FS"),
              file.path(keygenes_dir, "featureplots"))) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# 设置工作目录
setwd(resultdir)

# ======================================= 函数引入 =======================================
source(file.path(scriptdir, "Functions_scRNA.R"))
source(file.path(scriptdir, "Functions_Enrich.R"))

# ======================================= 加载数据 =======================================
# 读取 Seurat 对象
seurat_object <- readRDS(file.path(analysis_dir, "seurat_celltype.rds"))
setwd(keygenes_dir)

# ======================================= 基因集定义 =======================================
# 定义不同的基因集及其配置（可以方便地添加或修改）
gene_sets <- list(
  # 节律基因 (Circadian Rhythm)
  CR = list(
    genes = c("Arntl", "Clock", "Cry1", "Cry2", "Csnk1d", "Csnk1e", "Npas2", "Nr1d1", "Per1", "Per2", "Per3"),
    vln_filename = "Vlnplot_CR",
    feature_prefix = "CR",
    boxpair_dir = "boxpair_plots/CR",
    ncol = 4
  ),
  # 铁硫簇组装基因 (Fe-S Cluster)
  FS = list(
    genes = c("Abcb7", "AK157302", "Bola2", "Bola3", "Ciao1", "Ciao2a", "Ciao2b", "Ciao3", 
              "Ciapin1", "Fdx2", "Fxn", "Glrx3", "Glrx5", "Hscb", "Hspa9", "Iba57", "Isca1", 
              "Isca2", "Iscu", "Lyrm4", "Ndor1", "Ndufab1", "Ndufab1-ps", "Nfs1", "Nfu1", 
              "Nubp1", "Nubp2", "Nubpl", "Xdh"),
    vln_filename = "Vlnplot_FS",
    feature_prefix = "FS",
    boxpair_dir = "boxpair_plots/FS",
    ncol = 6
  )
)

# ======================================= 绘图配置 =======================================
# 绘图参数
plot_params <- list(
  group_colors = c(D = "#dd5633", L = "#5fb9d4"),
  group_var = "ftype",
  condition_var = "mfield"
)

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

# ======================================= 配对箱线图 =======================================
# 根据group变量分组
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

# ======================================= 基因表达图 =======================================
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

# 输出配置信息用于记录
config_info <- list(
  project = project_name,
  resolution = clustering_resolution,
  gene_sets = gene_sets,
  plot_params = plot_params,
  directories = list(
    base = base_dir,
    results = resultdir,
    project = project_dir,
    analysis = analysis_dir,
    keygenes = keygenes_dir
  )
)

saveRDS(config_info, file.path(keygenes_dir, "analysis_config.rds"))
