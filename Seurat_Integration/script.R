
##########----------------READY? START!--------------################

setwd("HCC_fetal/")

library(Seurat)
library(dplyr)
library(stringr)
library(scales)
library(ggplot2)
library(tidyr)

wholeatlas <- readRDS('/mnt/justinedata/HCC_fetal/data/wholeatlas.RData')


##########----------------formatting data--------------################

HCC <- ReadH5AD(file = "/home/rstudio/rawdata/Liver10X_73k.h5ad")
fetal4 <- ReadH5AD(file = "/home/rstudio/rawdata/fetal_8x_v1.h5ad")

head(HCC@meta.data)
unique(HCC@meta.data$orig.ident)
HCC@meta.data['tech']='HCC'
Idents(object=HCC) <- HCC@meta.data$tech
head(HCC@active.ident)
HCC@assays$RNA@counts <- HCC@assays$RNA@data
saveRDS(HCC, file='/home/rstudio/rawdata/Liver10X_73k.RData')

head(fetal4@meta.data)
unique(fetal4@meta.data$orig.ident)
fetal4@meta.data['tech']='fetal'
Idents(object=fetal4) <- fetal4@meta.data$tech
head(fetal4@active.ident)
fetal4@assays$RNA@counts <- fetal4@assays$RNA@data
saveRDS(fetal4, file='/home/rstudio/rawdata/fetal4.RData')

#####################################################################################


#############---------------------Integration-------------------################

merged <- merge(HCC, y = fetal4, add.cell.ids = c("HCC",'fetal'), project = 'Liver')
merged.list <- SplitObject(object=merged, split.by='tech')
merged.list <- merged.list[c('HCC','fetal')]

for (i in 1:length(merged.list)) {
  merged.list[[i]] <- NormalizeData(merged.list[[i]], verbose = FALSE)
  merged.list[[i]] <- FindVariableFeatures(merged.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

reference.list <- merged.list[c('HCC','fetal')]

reference.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

reference.integrated <- IntegrateData(anchorset = reference.anchors, dims = 1:30)

reference.integrated <- ScaleData(object = reference.integrated, verbose=FALSE)

reference.integrated <- RunPCA(object = reference.integrated, npcs = 100, verbose = FALSE, do.print=TRUE,
                               pcs.print = 1:5, genes.print=5)

ElbowPlot(reference.integrated, ndims = 80)

reference.integrated <- FindNeighbors(reference.integrated, dims = 1:50)

reference.integrated <- FindClusters(reference.integrated, resolution = 0.85)

reference.markers <- FindAllMarkers(reference.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

reference.integrated <- RunUMAP(reference.integrated, reduction = "pca", dims = 1:50)

saveRDS(wholeatlas, file = "/mnt/cellrangercount/HCC/HCC_fetal/data/wholedata.RData")


#####################################################################################


##############--------------------------PLOTS------------------------#################

DimPlot(reference.integrated, reduction = "umap", group.by = "tech", pt.size=0.8)
DimPlot(reference.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

FeaturePlot(reference.integrated, min.cutoff = 0, pt.size = 1,
            features = c('ALB',"ACTA2",'PLVAP','VWF','IL7R','CD8A','S100A8','CD163','MZB1','LYZ'))

FeaturePlot(reference.integrated, min.cutoff = 0, pt.size = 0.6,
            features = c('CD3D','ACTA2','VWF','MZB1','LYZ','KRT18'))


top10 <- reference.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
pl <- DoHeatmap(reference.integrated, features = top10$gene) + NoLegend()

pdf("figures/heatmap_pc50_cl47.pdf", width=30, height=30, onefile=T)
pl
dev.off()


DimPlot(wholeatlas, reduction = "umap", group.by = "PatientID", pt.size=0.8)

#####################################################################################




##########----------------Reformating columns--------------################

# --> when loading obj from h5ad, some information reformat to numerical 
# --> some double checking is advisable 

df <- wholeatlas@meta.data[c('tech','fetal.age','NormalvsTumor','PatientID')]

# reformat normal tumor 
# 1 --> Tumor
# 0 --> Normal
df['NormalvsTumor'] <- df$NormalvsTumor %>% 
  str_replace_all('1', "Tumor") %>% 
  str_replace_all('0', "Adj Normal")
unique(df$NormalvsTumor)


# reformat fetal age
# 0 --> Fw4
# 1 --> Fw6
# 2 --> Fw8
# 3 --> Fw21
df['fetal.age'] <- df$fetal.age %>% 
  str_replace('0', "Fw4") %>% 
  str_replace('1', "Fw6") %>% 
  str_replace('2', "Fw8") %>% 
  str_replace('3', "Fw21")

df['fetal.age'] <- df$fetal.age %>% 
  str_replace('Fw4', "Fw14") %>% 
  str_replace('Fw6', "Fw16") %>% 
  str_replace('Fw8', "Fw18")
unique(df$fetal.age)

# merging Normal Tumor and fetal age 
df1 <- df %>% mutate(mycol = coalesce(fetal.age, NormalvsTumor)) %>%
  select( mycol)
rownames(df1) <- rownames(df)
reference.integrated@meta.data['condition'] <-df1

# merging Normal Tumor and fetal
table(df$tech)

df$tech[df$tech %in% c("HCC")] <- NA
head(df)
df1 <- df %>% mutate(mycol = coalesce(tech, NormalvsTumor)) %>%
  select( mycol)
rownames(df1) <- rownames(df)
reference.integrated@meta.data['NTF'] <-df1
table(df1$mycol)

## rename patientID

wholeatlas@meta.data['PatientID'] <- paste0('P', wholeatlas$patientno)
library(car)
wholeatlas@meta.data['PatientID'] <-recode(wholeatlas$PatientID,"c('P0')=c('HN')")
wholeatlas@meta.data['PatientID'] <-recode(wholeatlas$PatientID,"c('PNA')=c('NA')")

# --> not sure why mutate does not work here :'(
if(FALSE) {
  df1 <- df %>% mutate(mycol = coalesce(PatientID, fetal.age)) %>%
    select( mycol)
  rownames(df1) <- rownames(df)
  reference.integrated@meta.data['PatientID'] <-df1
}

# --> easy way out paste columns tgt and remove NAs
wholeatlas@meta.data['fetal.age'] <- df$fetal.age
wholeatlas@meta.data['PatientID'] <- paste(wholeatlas$PatientID, wholeatlas$fetal.age, sep = "")
wholeatlas@meta.data['PatientID'] <- gsub("NA", "", wholeatlas$PatientID)

table(wholeatlas$PatientID)

# --> and TAADAAAAAA
DimPlot(wholeatlas, reduction = "umap", group.by = "PatientID", pt.size=0.8)

#####################################################################################


##########----------------PLOTS--------------################

reference.markers <- readRDS('/mnt/justinedata/HCC_fetal/data/reference.markers.RData')

top10 <- reference.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
DotPlot(reference.integrated, features = top10$gene[!duplicated(top10$gene)],
        col.min = -1, col.max = 2.5, dot.scale = 5) + 
  scale_colour_gradientn(colours = rev(viridis)) + RotatedAxis()

DimPlot(reference.integrated, reduction = "umap", group.by = "condition", pt.size=0.8,
        cols = c( "#1089eb", "#00ccb1", "#b52ec9", "#FF7F00", "#14690a", "#a62626"))

colorcode <- c('Adj Normal'= alpha("#1089eb", 0.3), 'fetal'= alpha("#00ccb1", 0.3), 'Tumor'= alpha("#a62626", 0.3))
DimPlot(reference.integrated, reduction = "umap", group.by = "NTF", pt.size=0.5,
        cols = colorcode)

FeaturePlot(reference.integrated, min.cutoff = 0, pt.size = 0.1, blend.threshold = 0.01, 
            cols = magma, features = c('ALB','KRT8','KRT19'), ncol = 2)

df <- reference.markers[,c('cluster', 'gene')]
rankgenes <- as.data.frame( df %>% 
                              group_by(cluster) %>%  # group by everything other than the value column. 
                              mutate(row_id=1:n()) %>% ungroup() %>%  # build group index
                              spread(key=cluster, value=gene) %>%    # spread
                              select(-row_id))  # drop the index
rankgenes[c(1:100),].to_csv(file='top100genes_39cl.csv')

genes = c('ALB', 'KRT19', 'PLVAP', 'ACTA2', 'IL7R', 'CD8A', 'GNLY', 'MZB1', 'CD14')
for (i in genes){
  pdf(paste0("figures/wholeatlas/", toString(i),".pdf"), width=7, height=5.7)
  print(FeaturePlot(reference.integrated, min.cutoff = 0, pt.size = 0.5, features = i, 
                    blend.threshold = 0.01, cols = magma))
  dev.off()
  
}

#####################################################################################


###########----------------COLORS--------------#################

magma = colorRampPalette(c('#050238','#59139e','#fa7916', '#fffdbd'), alpha = TRUE)(8)
viridis = colorRampPalette(c('#fde725','#35b778','#3e4a89','#430154'), alpha = TRUE)(8)

magma = c("#050238FF", "#280963FF", "#4D108FFF", "#863077FF", "#CC5B3CFF", "#FA8B2DFF", "#FCC475FF", "#FFFDBD1A")
magma = c("#0502381A", "#280963FF", "#4D108FFF", "#863077FF", "#CC5B3CFF", "#FA8B2DFF", "#FCC475FF", "#FFFDBDFF")

NTF_colors = c('#81007f', '#f5bb06', '#e6194b')

# Load the "scales" package
require(scales)
# Create vector with levels of object@ident
identities <- unique(mye$PatientID)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

################################################################


##########-----------Proportion barplot---------################

meta <- wholeatlas@meta.data[,colnames(wholeatlas@meta.data) %in% c('seurat_clusters','PatientID')]
head(meta)

df <- meta %>% group_by(seurat_clusters,PatientID) %>% summarise(n = n()) %>% spread(PatientID, n, fill = 0)
df <- as.data.frame(df)
rownames(df) <- df$seurat_clusters
df <- df[, -which(names(df) %in% c("seurat_clusters"))]
df <- df/colSums(df)*100
df <- df/rowSums(df)*100

par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(as.matrix(t(df)), horiz=TRUE, col = my_color_palette, legend = colnames(df),  args.legend = list(x = "topright", bty = "n", inset = c(-0.25, 0)), border = FALSE)


meta <- wholeatlas@meta.data[,colnames(wholeatlas@meta.data) %in% c('seurat_clusters','NTF')]
head(meta)

df <- meta %>% group_by(seurat_clusters,NTF) %>% summarise(n = n()) %>% spread(NTF, n, fill = 0)
df <- as.data.frame(df)
rownames(df) <- df$seurat_clusters
df <- df[, -which(names(df) %in% c("seurat_clusters"))]
df <- df/colSums(df)*100
df <- df/rowSums(df)*100

par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(as.matrix(t(df)), horiz=TRUE, col = NTF_colors , legend = colnames(df),  
        args.legend = list(x = "topright", bty = "n", inset = c(-0.25, 0)), border = FALSE)


################################################################

pdf(paste0("/mnt/justinedata/HCC_fetal/figures/wholeatlas/umapNTF_39.pdf"), width=9, height=5.7)
DimPlot(wholeatlas, reduction = "umap", group.by = 'NTF',cols = c('#81007f', '#f5bb06', '#e6194b'),  pt.size=0.5)
dev.off()

##########-----------rename cluster to immune non immune---------################
immune = c('36','0','17','6','7','27','3','16','38','37','18','13','30','23',
           '12','32','5','21','19','4','1','2')
nonimmune = c('14','11','29','25','24','35','33','22','26','31','8','15','34',
              '10','9','20',)

wholeatlas@meta.data['CD45'] <-wholeatlas@meta.data['seurat_clusters']

library(car)
wholeatlas@meta.data['CD45'] <-recode(wholeatlas$CD45,"c('36','0','17','6','7','27','3','16','38','37','18','13','30','23','12','32','5','21','19','4','1','2','28')=c('CD45pos')")
wholeatlas@meta.data['CD45'] <-recode(wholeatlas$CD45,"c('14','11','29','25','24','35','33','22','26','31','8','15','34','10','9','20')=c('CD45neg')")
wholeatlas@meta.data['NTF'] <-recode(wholeatlas$NTF,"c('fetal')=c('Fetal')")
wholeatlas$NTF <- factor(wholeatlas$NTF,
                       levels=c("Adj Normal","Tumor","Fetal" ))
wholeatlas$NTF[is.na(wholeatlas$NTF)] <- 'Adj Normal'

DimPlot(wholeatlas, reduction = "umap", group.by = "CD45", pt.size=0.8)

##########Proportion barplot

meta <- wholeatlas@meta.data[,colnames(wholeatlas@meta.data) %in% c('NTF','CD45')]
head(meta)

df <- meta %>% group_by(CD45,NTF) %>% summarise(n = n()) %>% spread(NTF, n, fill = 0)
df <- as.data.frame(df)
rownames(df) <- df$CD45
df <- df[, -which(names(df) %in% c("CD45"))]
df <- df/colSums(df)*100
df <- df/rowSums(df)*100

par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(as.matrix(t(df)), horiz=FALSE, col = c('#5000ff','#e6194B','#d6b8f5'), legend = colnames(df), beside=TRUE,  
        args.legend = list(x = "topright", bty = "n", inset = c(-0.4, 0)), border = FALSE)

##

df <- meta %>% group_by(NTF,CD45) %>% summarise(n = n()) %>% spread(CD45, n, fill = 0)
df <- as.data.frame(df)
rownames(df) <- df$NTF
df <- df[, -which(names(df) %in% c("NTF"))]
df <- df/colSums(df)*100
df <- df/rowSums(df)*100

par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(as.matrix(t(df)), horiz=FALSE, col = c('#010a43','#ffa41b'), legend = colnames(df), beside=TRUE,  
        args.legend = list(x = "topright", bty = "n", inset = c(-0.4, 0)), border = FALSE)

################################################





