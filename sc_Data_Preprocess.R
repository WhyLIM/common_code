
library(Seurat)
library(ggsc)
library(tidyverse)
library(ggpubr)
library(ggrepel)

# 清理环境
rm(list = ls())
gc()

# ======================================= 参数设置（可自定义） =======================================
# 基础目录设置 - 只需修改这个根目录，其他目录会自动生成
base_dir <- "~/scLC"
project_name <- "W_Penis"  # 项目名称，可以修改
sample_info_file <- "sample_ID.csv"  # 样本信息文件

# 样本选择条件
sample_conditions <- list(
  aging = "W",
  organ = "Penis"
)

# 质控和过滤参数
qc_params <- list(
  min_features = 300,
  max_features = 7500,
  max_counts = 50000,
  max_mito_percent = 15
)

# 降维参数
dim_reduction_params <- list(
  npcs = 20,
  dims_use = 1:20
)

# 聚类参数
clustering_resolutions <- c(0.3, 0.5, 0.8)

# 子目录自动生成
datadir <- file.path(base_dir, "data")
resultdir <- file.path(base_dir, "results")
scriptdir <- file.path(base_dir, "scripts")
project_dir <- file.path(resultdir, project_name)

# 自动创建项目目录
if (!dir.exists(project_dir)) {
  dir.create(project_dir, recursive = TRUE)
}

# 设置工作目录
setwd(resultdir)

# ======================================= 函数引入 =======================================
source(file.path(scriptdir, "Functions_scRNA.R"))

# ======================================= 样本加载 =======================================
# 读取样本信息
sample_ID <- read.csv(file.path(datadir, sample_info_file)) %>% 
  separate(abbr, into = c("mfield", "aging", "organ"), sep = "_", remove = FALSE)

# 根据条件过滤样本
filtered_samples <- sample_ID
for (cond_name in names(sample_conditions)) {
  filtered_samples <- filtered_samples %>% 
    filter(!!sym(cond_name) == sample_conditions[[cond_name]])
}

# 构建样本路径
sample_paths <- file.path(datadir, filtered_samples$production_ID)
names(sample_paths) <- filtered_samples$production_ID

# ======================================= 创建Seurat对象 =======================================
message("创建Seurat对象...")
seurat_list <- list()

for (name in names(sample_paths)) {
  path <- sample_paths[[name]]
  message(paste("处理样本:", name, "路径:", path))
  
  counts <- Read10X(data.dir = path, gene.column = 1)
  seurat_list[[name]] <- CreateSeuratObject(counts, project = name, min.cells = 3, min.features = 100)
}

seurat_object <- merge(x = seurat_list[[1]], y = seurat_list[-1])

# 清理内存
rm(seurat_list)
gc()

# ======================================= 质控和过滤 =======================================
message("执行质控和过滤...")
# 基本质控
seurat_object <- basic_qc(seurat_object, dir = project_dir)

# 构建过滤条件
filter_expr <- paste0(
  "nFeature_RNA > ", qc_params$min_features, 
  " & nFeature_RNA < ", qc_params$max_features,
  " & nCount_RNA < ", qc_params$max_counts,
  " & percent_mito < ", qc_params$max_mito_percent
)

# 应用过滤
seurat_object <- subset(seurat_object, subset = eval(parse(text = filter_expr)))

# 生成过滤后的质控图
p <- VlnPlot(seurat_object, group.by = "orig.ident", 
             features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb"), 
             ncol = 3, pt.size = 0)

width_calc <- length(unique(seurat_object$orig.ident))/2 + 10
ggsave(file.path(project_dir, "Vlnplot_mito_ribo_hb_filter.png"), p, width = width_calc, height = 10)
ggsave(file.path(project_dir, "Vlnplot_mito_ribo_hb_filter.pdf"), p, width = width_calc, height = 10)

# ======================================= 降维分析 =======================================
message("执行数据标准化和降维...")
# 标准化和降维
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object, verbose = FALSE)
seurat_object <- RunPCA(seurat_object, npcs = dim_reduction_params$npcs, verbose = FALSE)

# 批次效应校正
seurat_object <- IntegrateLayers(seurat_object, method = HarmonyIntegration, 
                                 orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
seurat_object <- JoinLayers(seurat_object)

# 可视化批次效应
p1 <- DimPlot(seurat_object, reduction = "pca", group.by = "orig.ident")
ggsave(file.path(project_dir, "batch_effect_pca.png"), p1, width = 7, height = 5)

p2 <- DimPlot(seurat_object, reduction = "harmony", group.by = "orig.ident")
ggsave(file.path(project_dir, "batch_effect_harmony.png"), p2, width = 7, height = 5)

# ======================================= 聚类和UMAP =======================================
message("执行聚类和UMAP降维...")
# 寻找邻居
seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = dim_reduction_params$dims_use)

# 多分辨率聚类
for (res in clustering_resolutions) {
  message(paste("执行聚类，分辨率:", res))
  seurat_object <- FindClusters(seurat_object, resolution = res)
}

# 运行UMAP
seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = dim_reduction_params$dims_use)

# ======================================= 添加元数据 =======================================
message("添加样本元数据...")
cmeta <- seurat_object@meta.data
# 将 cmeta 的行名保存为一个新列
cmeta <- cmeta %>% mutate(row_names = rownames(cmeta))
cmeta <- cmeta %>%
  left_join(filtered_samples, by = c("orig.ident" = "production_ID")) %>% 
  column_to_rownames("row_names")
seurat_object@meta.data <- cmeta

# 转换因子变量
factor_cols <- c("abbr", "mfield", "aging", "organ")
for (col in factor_cols) {
  if (col %in% colnames(seurat_object@meta.data)) {
    seurat_object@meta.data[[col]] <- as.factor(seurat_object@meta.data[[col]])
  }
}

# ======================================= 可视化和保存结果 =======================================
message("生成可视化结果并保存...")
# 为每个聚类分辨率生成UMAP图
for (res_col in grep("RNA_snn_res", colnames(seurat_object@meta.data), value = TRUE)) {
  res_val <- sub(".*res\\.", "", res_col)
  p <- DimPlot(seurat_object, reduction = "umap", group.by = res_col, label = TRUE) +
    ggtitle(paste("Resolution:", res_val))
  ggsave(file.path(project_dir, paste0("umap_res", res_val, ".png")), p, width = 10, height = 8)
}

# 按样本可视化
p <- DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident")
ggsave(file.path(project_dir, "umap_by_sample.png"), p, width = 10, height = 8)

# 按其他因素可视化
for (factor in c("mfield", "aging", "organ")) {
  if (factor %in% colnames(seurat_object@meta.data)) {
    p <- DimPlot(seurat_object, reduction = "umap", group.by = factor)
    ggsave(file.path(project_dir, paste0("umap_by_", factor, ".png")), p, width = 10, height = 8)
  }
}

# 保存Seurat对象
saveRDS(seurat_object, file.path(project_dir, paste0("seurat_", project_name, ".rds")))

# 保存配置信息用于记录
config_info <- list(
  project = project_name,
  sample_conditions = sample_conditions,
  qc_params = qc_params,
  dim_reduction_params = dim_reduction_params,
  clustering_resolutions = clustering_resolutions,
  directories = list(
    base = base_dir,
    data = datadir,
    results = resultdir,
    project = project_dir
  )
)

saveRDS(config_info, file.path(project_dir, "analysis_config.rds"))

message("分析完成！结果保存在: ", project_dir)