setwd('/mnt/justinedata/nichenetr/hepa_endo/')

library(nichenetr)
library(Seurat)
library(tidyverse)
library(ggplot2)

wholeatlas <- readRDS("/mnt/justinedata/HCC_fetal/data/reference.integrated.RData")

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5]

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr)


wholeatlas@meta.data['DC1'] <- apply(FetchData(wholeatlas, c('CADM1', 'XCR1', 'CLEC9A', 'CD74')),1, median)
wholeatlas@meta.data['DC2'] <- apply(FetchData(wholeatlas, c('HLA-DRA', 'CD1C', 'FCER1A', 'HLA-DPB1', 'CLEC10A')),1, median)
wholeatlas@meta.data['pDC'] <- apply(FetchData(wholeatlas, c('IRF7', 'LILRA4', 'IL3RA', 'IGKC', 'BCL11A', 'GZMB')),1, median)

magma = c("#0502384D", "#280963FF", "#4D108FFF", "#863077FF", "#CC5B3CFF", "#FA8B2DFF", "#FCC475FF", "#FFFDBDFF")
FeaturePlot(wholeatlas, min.cutoff = 0, pt.size = 0.8, blend.threshold = 0.01, cols = magma,
            features = c('DLL4','IL6'))


df <- as.data.frame(wholeatlas@meta.data$seurat_clusters)
colnames(df) = "CellType"

library(car)
df$CellType <-recode(df$CellType,"c(14,11,29,24,25,35)=c('Hepatocytes')")
df$CellType <-recode(df$CellType,"c(28)=c('Bcells')")
df$CellType <-recode(df$CellType,"c(37,18,23,13,30)=c('NA')")
df$CellType <-recode(df$CellType,"c(32,5,19,21,4,1,2)=c('Erythroids')")
df$CellType <-recode(df$CellType,"c(22,26,31)=c('Fibroblast')")
df$CellType <-recode(df$CellType,"c(8,15,34,10,9,20)=c('Endothelial')")
df$CellType <-recode(df$CellType,"c(36,17,0,6,7,27,3,16)=c('Immune')")
df$CellType <-recode(df$CellType,"c(38,33)=c('Others')")
df$CellType <-recode(df$CellType,"c(12)=c('Fetal B')")

unique(df$CellType)

wholeatlas@meta.data["CellType"] <- df$CellType

DimPlot(mye, reduction = "umap",  label = TRUE, pt.size = 0.5) + NoLegend()#group.by = 'CellType',

Idents(object=wholeatlas) <- wholeatlas@meta.data$CellType

head(wholeatlas@meta.data)

wholeatlas@meta.data$CellType %>% table() 

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = wholeatlas, 
  receiver = "Endothelial", 
  condition_colname = "NTF", condition_oi = "Tumor", condition_reference = "Adj Normal", 
  sender = c("Hepatocytes"), 
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  organism = "human")

nichenet_output$ligand_activity_target_heatmap


nichenet_output$ligand_activities
nichenet_output$top_ligands

DotPlot(wholeatlas, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

nichenet_output$ligand_target_heatmap

nichenet_output$ligand_activity_target_heatmap

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = wholeatlas, 
  receiver = c("Myeloids"), 
  condition_colname = "NTF", condition_oi = "Tumor", condition_reference = "Adj Normal", 
  sender = c("Endothelial"), 
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  organism = "human")

nichenet_output$ligand_activity_target_heatmap

wholeatlas@


mye<- readRDS('/mnt/justinedata/HCC_fetal/data/mye.RData')

df <- as.data.frame(mye@meta.data$seurat_clusters)
colnames(df) = "CellType"

library(car)
df$CellType <-recode(df$CellType,"c(1)=c('FTAM')")
df$CellType <-recode(df$CellType,"c(7)=c('TAM2')")
df$CellType <-recode(df$CellType,"c(8,5,30)=c('FLM')")
df$CellType <-recode(df$CellType,"c(4)=c('TAM')")
df$CellType <-recode(df$CellType,"c(6,2,10)=c('Monocytes')")
df$CellType <-recode(df$CellType,"c(0,3,9,11)=c('DC')")

df$CellType <-recode(df$CellType,"c(1,4,7)=c('TAM')")
df$CellType <-recode(df$CellType,"c(8,5,30)=c('FLM')")
df$CellType <-recode(df$CellType,"c(6,2,10)=c('Monocytes')")
df$CellType <-recode(df$CellType,"c(0,3,9,11)=c('DC')")

df$CellType <-recode(df$CellType,"c(1,4,7,8,5,30)=c('TAM')")
df$CellType <-recode(df$CellType,"c(6,2,10)=c('Monocytes')")
df$CellType <-recode(df$CellType,"c(0,3,9,11)=c('DC')")


unique(df$CellType)

mye@meta.data["CellType"] <- df$CellType
DimPlot(mye, reduction = "umap", group.by='CellType', label = TRUE, pt.size = 0.5) + NoLegend()

df <- as.data.frame(wholeatlas@meta.data$seurat_clusters)
colnames(df) = "CellType"

library(car)
df$CellType <-recode(df$CellType,"c(14,11,29,24,25,35)=c('Hepatocytes')")
df$CellType <-recode(df$CellType,"c(28)=c('Bcells')")
df$CellType <-recode(df$CellType,"c(37,18,23,13,30)=c('NA')")
df$CellType <-recode(df$CellType,"c(32,5,19,21,4,1,2)=c('Erythroids')")
df$CellType <-recode(df$CellType,"c(22,26,31)=c('Fibroblast')")
df$CellType <-recode(df$CellType,"c(8,15,34,10,9,20)=c('Endothelial')")
df$CellType <-recode(df$CellType,"c(36,17,0,6,7,27,3,16)=c('Immune')")
df$CellType <-recode(df$CellType,"c(38,33)=c('Others')")
df$CellType <-recode(df$CellType,"c(12)=c('Fetal B')")

unique(df$CellType)

wholeatlas@meta.data["CellType"] <- df$CellType

DimPlot(wholeatlas, reduction = "umap", group.by = 'CellType', label = TRUE, pt.size = 0.5) + NoLegend()

myect <- mye@meta.data["CellType"]
wact <- wholeatlas@meta.data["CellType"]

df <- cbind(wact,myect[, "CellType"][match(rownames(wact), rownames(myect))])
colnames(df) <- c('wa','mye')
View(df)

df[df == 'NA'] <- NA

df['comb'] <- paste(df$wa,df$mye, sep = "")

df['comb'] <- gsub("NA", "", df$comb)

wholeatlas@meta.data["CT_mye"] <- df$comb

unique(df$comb)

DimPlot(wholeatlas, reduction = "umap", group.by='CT_mye', label = TRUE, pt.size = 0.5) + NoLegend()

Idents(object=wholeatlas) <- wholeatlas@meta.data$CT_mye



nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = wholeatlas, 
  receiver = c("TAM"), 
  condition_colname = "NTF", condition_oi = "Tumor", condition_reference = "Adj Normal", 
  sender = c("Endothelial"), 
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  organism = "human")

nichenet_output$ligand_activity_target_heatmap


plotnichenet_output = nichenet_seuratobj_cluster_de(
  seurat_obj = seuratObj, 
  receiver_reference = "CD8 T_SS", receiver_affected = "CD8 T_LCMV", 
  sender = c("DC_LCMV","Mono_LCMV"), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")


endo <- readRDS('/mnt/justinedata/HCC_fetal/data/endo.RData')

df <- FetchData(endo, c('DLL4', 'PLVAP'))
df['DLL4xPLVAP'] <- df$DLL4 > 1 & df$PLVAP > 1

endo@meta.data['DLL4xPLVAP'] <- df['DLL4xPLVAP']

DimPlot(endo, reduction = "umap", group.by = 'DLL4xPLVAP', pt.size = 0.5)#group.by = 'CellType',

VlnPlot(endo, "PLVAP", group.by = "NTF", pt.size = 0) + NoLegend()

FeaturePlot(endo,ncol = 1,min.cutoff = 0, pt.size = 0.8, blend.threshold = 0.01, cols = magma, features = c('DLL4','PLVAP'))

