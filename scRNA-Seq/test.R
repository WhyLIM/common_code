
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

source("~/Projects/common_code/scRNA-Seq/Functions_Enrich.R")

result <- plot_pathway(go_df = go_df, palette = 3, alpha = 0.5, 
                       top_n_per_category = 3, x_axis_column = "p.adjust")
result <- plot_pathway(go_df = go_df, kegg_df = kegg_df, palette = 3, alpha = 0.5, 
                       top_n_per_category = 3, x_axis_column = "FoldEnrichment")
result <- plot_pathway(go_df = go_df, kegg_df = kegg_df, palette = 3, alpha = 0.5, 
                       top_n_per_category = 3, x_axis_column = "GeneRatio")
result <- plot_pathway(go_df = go_df, kegg_df = kegg_df, palette = 3, alpha = 0.5, 
                       top_n_per_category = 3, x_axis_column = "zScore")
plot_pathway(kegg_df = kegg_df, top_n_per_category = 10)


