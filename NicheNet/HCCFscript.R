
HCCF <- ReadH5AD(file = "/mnt/justinedata/rawdata/adata_final_PC55.h5ad")

DimPlot(HCCF, reduction = "umap", group.by = "louvain", label = TRUE, pt.size = 2) + NoLegend()


df <- as.data.frame(HCCF@meta.data$louvain)
colnames(df) = "CellType"

library(car)
df$CellType <-recode(df$CellType,"c(5,22,27,25)=c('Myeloid')")
df$CellType <-recode(df$CellType,"c(21,9)=c('Bcell')")
df$CellType <-recode(df$CellType,"c(11,13,18,20,24,28)=c('Hepatocytes')")
df$CellType <-recode(df$CellType,"c(29)=c('Bipotent')")
df$CellType <-recode(df$CellType,"c(19,23)=c('Fibroblast')")
df$CellType <-recode(df$CellType,"c(2,8,10)=c('Endothelial')")
df$CellType <-recode(df$CellType,"c(0,16,17)=c('CD4')")
df$CellType <-recode(df$CellType,"c(3,12)=c('NK')")
df$CellType <-recode(df$CellType,"c(1,14,15)=c('CD8')")
df$CellType <-recode(df$CellType,"c(30)=c('Mast')")
df$CellType <-recode(df$CellType,"c(4,6,7)=c('Fetal Erythroid')")
df$CellType <-recode(df$CellType,"c(26)=c('Megakaryocytes')")
unique(df$CellType)

HCCF@meta.data["CellType"] <- df$CellType
DimPlot(HCCF, reduction = "umap", group.by = "CellType", label = TRUE, pt.size = 2) + NoLegend()

Idents(HCCF) <- HCCF@meta.data$CellType



df <- as.data.frame(HCCF@meta.data$NormalvsTumor)
colnames(df) = "NTF"

df$NTF <-recode(df$NTF,"c(0,1)=c('Fetal')")
df$NTF <-recode(df$NTF,"c(2)=c('Adj Normal')")
df$NTF <-recode(df$NTF,"c(3)=c('Tumor')")
unique(df$NTF)

HCCF@meta.data["NTF"] <- df$NTF
DimPlot(HCCF, reduction = "umap", group.by = "NTF", pt.size = 2)




mye <- ReadH5AD(file = "/mnt/justinedata/rawdata/mye_wo_b.h5ad")

DimPlot(mye, reduction = "umap", group.by = "louvain", label = TRUE, pt.size = 2) + NoLegend()

df <- as.data.frame(mye@meta.data$louvain)
colnames(df) = "CellType"

library(car)
df$CellType <-recode(df$CellType,"c(0)=c('TAM2')")
df$CellType <-recode(df$CellType,"c(1,10)=c('DC2')")
df$CellType <-recode(df$CellType,"c(2)=c('Monocytes')")
df$CellType <-recode(df$CellType,"c(3,9)=c('TAM1')")
df$CellType <-recode(df$CellType,"c(4)=c('FMono')")
df$CellType <-recode(df$CellType,"c(5)=c('FDC2')")
df$CellType <-recode(df$CellType,"c(6)=c('FLM')")
df$CellType <-recode(df$CellType,"c(7)=c('DC1')")
df$CellType <-recode(df$CellType,"c(8)=c('TAM3')")
df$CellType <-recode(df$CellType,"c(11)=c('pDC')")
df$CellType <-recode(df$CellType,"c(12)=c('Megakaryocytes')")
df$CellType <-recode(df$CellType,"c(13)=c('FCMP')")
unique(df$CellType)

mye@meta.data["CellType"] <- df$CellType
DimPlot(mye, reduction = "umap", group.by = "CellType", label = TRUE, pt.size = 2) + NoLegend()

myect <- mye@meta.data["CellType"]
wact <- HCCF@meta.data["CellType"]

df <- cbind(wact,myect[, "CellType"][match(rownames(wact), rownames(myect))])
colnames(df) <- c('wa','mye')
View(df)
df['comb'] <- paste(df$wa,df$mye, sep = "")
df['comb'] <- gsub("NA", "", df$comb)
df['comb'] <- gsub("Myeloid", "", df$comb)
df$comb[df$comb == ""] <-'pDC' 

HCCF@meta.data["CT_mye"] <- df$comb
DimPlot(HCCF, reduction = "umap", group.by='CT_mye', label = TRUE, pt.size = 0.5) + NoLegend()

Idents(HCCF) <- HCCF@meta.data$CT_mye





nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = HCCF, 
  receiver = "TAM1", 
  condition_colname = "NTF", condition_oi = "Tumor", condition_reference = "Adj Normal", 
  sender = c("Endothelial"), 
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  organism = "human")

nichenet_output$ligand_activity_target_heatmap



FeaturePlot(HCCF, min.cutoff = 0, pt.size = 0.5, ncol=2, cols = magma,
            features = c('DLL4','HES1','PLVAP','FOLR2'))










































