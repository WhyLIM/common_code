
library(Seurat)
library(ggsc)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggtangle)
library(ggalluvial)
library(openxlsx)
library(rlang)
library(dplyr)

# 清理环境
rm(list = ls())
gc()

# ======================================= 参数设置 =======================================
# 基础目录设置
base_dir <- "~/scLC"
project_name <- "W_Penis"  # 项目名称

seurat_obj <- "seurat_celltype.rds"
clustering_resolution <- "0.5"  # 聚类分辨率
group_var <- "ftype"  # 分组变量，如设置为 celltype 则按细胞类型来分别进行差异分析
condition_var <- "mfield"  # 条件变量
ident_1 <- "L"  # 实验组
ident_2 <- "D"  # 对照组
deg_pval_cutoff <- 0.05  # 差异基因 P 值阈值
deg_top_n <- 200  # 富集分析时每组取前多少个基因
deg_celltype_top_n <- 10 # 每个细胞类型提取前多少个上下调基因（用于 ppt）
kegg_top_n <- 10  # KEGG 网络图显示前多少个通路
go_top_n <- 3  # GO 网络图每个类别显示前多少条目

# 自动生成子目录
datadir <- file.path(base_dir, "data")
resultdir <- file.path(base_dir, "results")
scriptdir <- file.path(base_dir, "scripts")
project_dir <- file.path(resultdir, project_name)
analysis_dir <- file.path(project_dir, paste0("check-markers-", clustering_resolution))
dea_dir <- file.path(analysis_dir, "DEA")
enrich_dir <- file.path(dea_dir, "Enrich")
kegg_dir <- file.path(enrich_dir, "KEGG")
go_dir <- file.path(enrich_dir, "GO")
keggnet_dir <- file.path(enrich_dir, "KEGG_Net")
gonet_dir <- file.path(enrich_dir, "GO_Net")
volcano_dir <- file.path(dea_dir, "volcano_plots")

# 创建所需目录
for (dir in c(dea_dir, enrich_dir, kegg_dir, go_dir, keggnet_dir, gonet_dir, volcano_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}


# ================================= 以下部分都不需要修改 =================================
# ======================================= 函数引入 =======================================
source(file.path(scriptdir, "Functions_DEA.R"))
source(file.path(scriptdir, "Functions_scRNA.R"))
source(file.path(scriptdir, "Functions_Enrich.R"))

# ======================================= 辅助函数 =======================================
# 创建带有多个工作表的 Excel 文件
create_excel_with_sheets <- function(data_list, file_name) {
  wb <- createWorkbook()
  for (sheet_name in names(data_list)) {
    safe_name <- substr(gsub("\\[|\\]|/|\\\\|\\?|\\*|:", "_", sheet_name), 1, 31)  # Excel 限制工作表名长度
    addWorksheet(wb, sheetName = safe_name)
    writeData(wb, sheet = safe_name, x = data_list[[sheet_name]])
  }
  saveWorkbook(wb, file = file_name, overwrite = TRUE)
  return(file_name)
}

# 过滤差异基因用于富集分析
filter_degs_for_enrichment <- function(degs, n_per_direction = 200) {
  lapply(degs, function(df) {
    # 按绝对 log2FC 值排序
    df <- df %>% 
      mutate(abs_log2FC = abs(avg_log2FC)) %>%
      arrange(desc(abs_log2FC))
    
    if (nrow(df) > 2 * n_per_direction) {
      # 分别提取上下调的前 n 个基因
      df <- df %>% 
        group_by(regulated) %>% 
        slice_head(n = n_per_direction) %>% 
        ungroup()
    }
    return(df)
  })
}

# 提取每个细胞类型中上下调的前 n 个基因
extract_top_degs <- function(degs, n = 10) {
  result <- lapply(names(degs), function(celltype) {
    df <- degs[[celltype]] %>%
      mutate(abs_log2FC = abs(avg_log2FC)) %>%
      arrange(desc(abs_log2FC))
    
    up <- df %>% filter(regulated == "Up") %>% slice_head(n = n)
    down <- df %>% filter(regulated == "Down") %>% slice_head(n = n)
  }) %>% bind_rows()
  
  return(result)
}

# KEGG 富集结果转换为网络图格式
format_kegg_for_network <- function(kegg_list, n) {
  result_list <- list()
  for (item_name in names(kegg_list)) {
    item <- kegg_list[[item_name]]
    top_kegg <- item$top_kegg
    top_n <- head(top_kegg, n)
    
    inner_list <- list()
    for (i in 1:nrow(top_n)) {
      description <- gsub(" - Mus musculus \\(house mouse\\)$", "", top_n$Description[i])
      geneID <- unlist(strsplit(top_n$geneID[i], "/"))
      inner_list[[description]] <- geneID
    }
    result_list[[item_name]] <- inner_list
  }
  return(result_list)
}

# GO 富集结果转换为网络图格式
format_go_for_network <- function(go_list, n) {
  result_list <- list()
  for (item_name in names(go_list)) {
    item <- go_list[[item_name]]
    top_go <- item$top_go
    
    # 按 ONTOLOGY 分组取前 n 行
    top_n <- top_go %>%
      group_by(ONTOLOGY) %>%
      slice_head(n = n) %>%
      ungroup()
    
    inner_list <- list()
    for (i in 1:nrow(top_n)) {
      description <- top_n$Description[i]
      geneID <- unlist(strsplit(top_n$geneID[i], "/"))
      inner_list[[description]] <- geneID
    }
    result_list[[item_name]] <- inner_list
  }
  return(result_list)
}

# 绘制保存网络图
save_network_plot <- function(cnet_data, prefix, title_prefix, dimensions = c(9, 8), show_category = 5, dir_path) {
  for (n in names(cnet_data)) {
    safe_name <- gsub("[/ ]", "_", n)
    filename_png <- file.path(dir_path, paste0(prefix, "_", safe_name, ".png"))
    filename_pdf <- file.path(dir_path, paste0(prefix, "_", safe_name, ".pdf"))
    
    p <- cnetplot(cnet_data[[n]], 
                  color_category = 'firebrick', 
                  color_item = 'steelblue', 
                  node_label = "gene", 
                  showCategory = show_category) + 
      geom_cnet_label(node_label = 'category', color='firebrick', size = 4) + 
      labs(title = paste(title_prefix, n)) + 
      theme(plot.margin = margin(t = 10, b = 10, l = 20, r = 20))
    
    ggsave(filename_png, plot = p, width = dimensions[1], height = dimensions[2])
    ggsave(filename_pdf, plot = p, width = dimensions[1], height = dimensions[2])
  }
}

# ======================================= 数据加载 =======================================
seurat_object <- readRDS(file.path(analysis_dir, seurat_obj))

# ===================================== 差异表达分析 =====================================
setwd(dea_dir)
# 1. 按细胞类型拆分 Seurat 对象
Idents(seurat_object) <- group_var
celltypes <- levels(seurat_object@meta.data[[group_var]])
seurat_cells <- list()

for (ct in celltypes) {
  seurat_cells[[ct]] <- subset(seurat_object, subset = !!sym(group_var) == ct)
}

# 2. 计算每个细胞类型中条件间的差异基因
message("计算差异表达基因...")
deg_list <- list()

for (celltype_name in names(seurat_cells)) {
  message(paste("分析细胞类型:", celltype_name))
  subset_obj <- seurat_cells[[celltype_name]]
  
  # 设置分组标识为条件变量
  Idents(subset_obj) <- condition_var
  
  # 计算差异基因
  markers <- FindMarkers(
    subset_obj,
    ident.1 = ident_1,   # 实验组
    ident.2 = ident_2,   # 对照组
    group.by = condition_var
  )
  
  # 标记上下调
  markers$regulated <- ifelse(markers$avg_log2FC > 0, "Up", "Down")
  
  # 存入列表
  deg_list[[celltype_name]] <- markers
}
# 保存差异基因结果
saveRDS(deg_list, "deg_list.rds")
# deg_list <- readRDS("deg_list.rds")

# 3. 合并所有细胞类型的差异基因为一个数据框
combined_degs <- do.call(rbind, lapply(names(deg_list), function(ct) {
  df <- deg_list[[ct]] %>% rownames_to_column(var = "SYMBOL")
  cbind(celltype = ct, df)
}))
# 保存合并的差异基因结果
write.csv(combined_degs, paste0("DEGs_celltype_", project_name, ".csv"), row.names = FALSE)
combined_degs <- read.csv(paste0("DEGs_celltype_", project_name, ".csv"))

# 4. 筛选显著差异基因
sig_markers <- combined_degs %>% filter(p_val_adj < deg_pval_cutoff)

# 5. 提取所有显著的上下调基因并按细胞类型保存为 Excel
up_genes <- sig_markers %>%
  filter(avg_log2FC > 0) %>%
  arrange(celltype, desc(abs(avg_log2FC))) %>%
  split(.$celltype)
down_genes <- sig_markers %>%
  filter(avg_log2FC < 0) %>%
  arrange(celltype, desc(abs(avg_log2FC))) %>%
  split(.$celltype)
create_excel_with_sheets(up_genes, "DEGs_celltype_up.xlsx")
create_excel_with_sheets(down_genes, "DEGs_celltype_down.xlsx")

# 6. 提取每个细胞类型中前 n 个上下调基因并保存
top_degs_df <- extract_top_degs(
  sig_markers %>% split(.$celltype),
  n = deg_celltype_top_n
)
write.csv(top_degs_df, paste0("DEGs_celltype_top", deg_celltype_top_n, ".csv"), row.names = FALSE)

# 7. 绘制火山图
message("绘制火山图...")
for (celltype_name in names(deg_list)) {
  markers <- deg_list[[celltype_name]]
  
  # 筛选显著基因标签
  labels <- markers %>%
    filter(abs(avg_log2FC) > 1 & p_val_adj < deg_pval_cutoff) %>% 
    rownames()
  
  # 绘制火山图
  p <- volcanoPlot(
    markers,
    logFC_col = "avg_log2FC",
    pvalue_col = "p_val_adj",
    labels = labels, 
    title = paste("Volcano plot of", gsub("[/ ]", "_", celltype_name))
  )
  
  # 保存图片
  filename <- paste0(volcano_dir, "/volcano_", project_name, "_", gsub("[/ ]", "_", celltype_name))
  ggsave(paste0(filename, ".png"), plot = p, width = 10, height = 8)
  ggsave(paste0(filename, ".pdf"), plot = p, width = 10, height = 8)
}

# ======================================= 富集分析 =======================================
setwd(enrich_dir)

# 1. 筛选用于富集分析的差异基因
message("准备富集分析...")
# 获取按细胞类型分组的差异基因
deg_celltype <- sig_markers %>% split(.$celltype)
# 按上下调和 logFC 过滤基因
degs_for_enrichment <- filter_degs_for_enrichment(deg_celltype, n_per_direction = deg_top_n)

# 2. KEGG 富集分析
message("执行 KEGG 富集分析...")
kegg_results <- list()

for (ct in names(degs_for_enrichment)) {
  DEG <- degs_for_enrichment[[ct]]
  
  # 上调基因 KEGG 富集
  DEG_up <- DEG %>% filter(regulated == "Up")
  safe_name <- gsub("[/ ]", "_", paste0(ct, "_up"))
  kegg_results[[paste0(ct, "_up")]] <- KEGGEnrichment(
    DEG_up, suffix = safe_name, species = "mouse", outputDir = kegg_dir
  )
  
  # 下调基因 KEGG 富集
  DEG_down <- DEG %>% filter(regulated == "Down")
  safe_name <- gsub("[/ ]", "_", paste0(ct, "_down"))
  kegg_results[[paste0(ct, "_down")]] <- KEGGEnrichment(
    DEG_down, suffix = safe_name, species = "mouse", outputDir = kegg_dir
  )
}

# 保存 KEGG 富集结果
saveRDS(kegg_results, "kegg_enrichment_results.rds")
# kegg_results <- readRDS("kegg_enrichment_results.rds")

# 3. GO 富集分析
message("执行 GO 富集分析...")
go_results <- list()

for (ct in names(degs_for_enrichment)) {
  DEG <- degs_for_enrichment[[ct]]
  
  # 上调基因 GO 富集
  DEG_up <- DEG %>% filter(regulated == "Up")
  safe_name <- gsub("[/ ]", "_", paste0(ct, "_up"))
  go_results[[paste0(ct, "_up")]] <- GOEnrichment(
    DEG_up, suffix = safe_name, species = "mouse", outputDir = go_dir
  )
  
  # 下调基因 GO 富集
  DEG_down <- DEG %>% filter(regulated == "Down")
  safe_name <- gsub("[/ ]", "_", paste0(ct, "_down"))
  go_results[[paste0(ct, "_down")]] <- GOEnrichment(
    DEG_down, suffix = safe_name, species = "mouse", outputDir = go_dir
  )
}

# 保存 GO 富集结果
saveRDS(go_results, "go_enrichment_results.rds")
# go_results <- readRDS("go_enrichment_results.rds")

# 4. 生成 KEGG 网络图
message("生成 KEGG 网络图...")
cnet_kegg <- format_kegg_for_network(kegg_results, kegg_top_n)
save_network_plot(
  cnet_kegg, 
  prefix = paste0("KEGGNet_", kegg_top_n), 
  title_prefix = "KEGG Network of",
  show_category = kegg_top_n,
  dir_path = keggnet_dir
)

# 5. 生成 GO 网络图
message("生成 GO 网络图...")
cnet_go <- format_go_for_network(go_results, go_top_n)
save_network_plot(
  cnet_go, 
  prefix = paste0("GONet_", go_top_n), 
  title_prefix = "GO Network of",
  show_category = go_top_n * 3,
  dir_path = gonet_dir
)
