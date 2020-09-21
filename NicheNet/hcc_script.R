hcc <- readRDS("/mnt/justinedata/rawdata/Liver10X_73k.RData") 

library(nichenetr)
library(Seurat)
library(tidyverse)
library(ggplot2)


DimPlot(hcc, reduction = "umap", group.by = "louvain", label = TRUE, pt.size = 2) + NoLegend()

df <- as.data.frame(hcc@meta.data$louvain)
colnames(df) = "CellType"


library(car)
df$CellType <-recode(df$CellType,"c('16','11','21','20','14','25','7','24')=c('Hepatocytes')")
df$CellType <-recode(df$CellType,"c(24)=c('Bipotent')")
df$CellType <-recode(df$CellType,"c(18)=c('Bcells')")
df$CellType <-recode(df$CellType,"c(13)=c('Fibroblast')")
df$CellType <-recode(df$CellType,"c(27,2,6,17,26)=c('Endothelial')")
df$CellType <-recode(df$CellType,"c(36,17,0,6,7,27,3,16)=c('Immune')")
df$CellType <-recode(df$CellType,"c(0,5,23)=c('CD4')")
df$CellType <-recode(df$CellType,"c(1,9,15,22)=c('NK')")
df$CellType <-recode(df$CellType,"c(3,4)=c('CD8')")
df$CellType <-recode(df$CellType,"c(28)=c('Mast')")
df$CellType <-recode(df$CellType,"c(12)=c('Treg')")
df$CellType <-recode(df$CellType,"c(8,10,19)=c('Myeloid')")
unique(df$CellType)

hcc@meta.data["CellType"] <- df$CellType
DimPlot(hcc, reduction = "umap", group.by = "CellType", label = TRUE, pt.size = 2) + NoLegend()

head(hcc@meta.data)

Idents(hcc) <- hcc@meta.data$CellType

head(hcc@meta.data)
df <- hcc@meta.data[c('NormalvsTumor','tech')]
df['NormalvsTumor'] <- df$NormalvsTumor %>% 
  str_replace_all('1', "Tumor") %>% 
  str_replace_all('0', "Adj Normal")
unique(df$NormalvsTumor)
hcc@meta.data["NormalvsTumor"] <- df$NormalvsTumor



nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = hcc, 
  receiver = "Myeloid", 
  condition_colname = "NormalvsTumor", condition_oi = "Tumor", condition_reference = "Adj Normal", 
  sender = c("Endothelial"), 
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  organism = "human")

nichenet_output$ligand_activity_target_heatmap
