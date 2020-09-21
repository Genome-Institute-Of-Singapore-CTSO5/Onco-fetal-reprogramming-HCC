
##########----------------READY? START!--------------################

mye <- readRDS('/mnt/cellrangercount/HCC/HCC_fetal/data/mye.RData')

mye <- subset(reference.integrated, idents = c(37,18,23,13,30))
mye

DimPlot(mye, reduction = "umap", label = TRUE, pt.size = 2) + NoLegend()


##########----------------Usual Practice--------------################

mye <- RunPCA(mye, features = VariableFeatures(object = mye), npcs = 70)

ElbowPlot(mye, ndims = 50)

mye <- FindNeighbors(mye, dims= 1:15)
mye <- RunUMAP(mye, dims = 1:15, min.dist=0.01, spread=0.3)

mye <- FindClusters(mye, resolution = 0.6)
DimPlot(mye, reduction = "umap", pt.size = 0.8, 
        label.size = 6, group.by = c('condition'))+ NoLegend()


saveRDS(mye, file = '/mnt/justinedata/HCC_fetal/data/mye.RData')

#####################################################################################


#############---------------------PLOTS-------------------################

cond = c('seurat_clusters', 'NTF', 'PatientID')

DimPlot(mye, reduction = "umap", pt.size=1, label=TRUE) + NoLegend()

DimPlot(mye, reduction = "umap", group.by = "NTF", pt.size=1, cols = c("#1089eb", "#00ccb1","#a62626"))

mye@meta.data['DC1'] <- apply(FetchData(mye, c('CADM1', 'XCR1', 'CLEC9A', 'CD74')),1, median)
mye@meta.data['DC2'] <- apply(FetchData(mye, c('HLA-DRA', 'CD1C', 'FCER1A', 'HLA-DPB1', 'CLEC10A')),1, median)
mye@meta.data['pDC'] <- apply(FetchData(mye, c('IRF7', 'LILRA4', 'IL3RA', 'IGKC', 'BCL11A', 'GZMB')),1, median)


FeaturePlot(mye, min.cutoff = 0, pt.size = 0.1, blend.threshold = 0.01, cols = magma, ncol = 3,
            features = c('S100A8','FOLR2','SPP1','MT1G','DC1','DC2','pDC'))

FeaturePlot(mye, min.cutoff = 0, pt.size = 0., blend.threshold = 0.01, cols = magma, ncol = 3,
            features = c('FOLR2','HES1'))

FeaturePlot(mye, min.cutoff = 0, pt.size = 0.5, features = i, 
            blend.threshold = 0.01, cols = magma)

VlnPlot(mye, i, group.by = "seurat_clusters", pt.size = 0) + NoLegend())


#####################################################################################


#############---------------------Renaming columns-------------------################

fetalage <-df[df$fetal.age %in% c("Fw14","Fw16","Fw18","Fw21"),'fetal.age',drop = FALSE]
NT <- df[df$NormalvsTumor %in% c("Tumor","Adj Normal"),'NormalvsTumor',drop = FALSE]
disease <- df[df$orig.ident %in% c("Normal","Cirrhotic"),'orig.ident',drop = FALSE]

df1 <- df %>% mutate(mycol = coalesce(fetal.age, NormalvsTumor, orig.ident)) %>%
  select( mycol)
rownames(df1) <- rownames(df)
unique(df1$mycol)
mye@meta.data['condition'] <-df1

mye@meta.data$condition[mye@meta.data$condition %in% c("Fw14","Fw16","Fw18","Fw21")] <- "Fetal"

mye <- mye[,!mye@meta.data$condition == 'Normal']

DimPlot(mye, reduction = "umap", group.by = "condition", pt.size=1)

##################################################################################


#############---------------------DE Genes-------------------################

# all against all
mye.markers <- FindAllMarkers(mye, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# differential genes between cluster 
cluster7.markers <- FindMarkers(mye, ident.1 = 7, ident.2 = 4, min.pct = 0.25)
head(cluster7.markers, n = 5)

#saving top100genes percluster 
df <- mye.markers[,c('cluster', 'gene')]
rankgenes <- as.data.frame( df %>% 
                              group_by(cluster) %>%  # group by everything other than the value column. 
                              mutate(row_id=1:n()) %>% ungroup() %>%  # build group index
                              spread(key=cluster, value=gene) %>%    # spread
                              select(-row_id))  # drop the index
write.csv(rankgenes[c(1:100),],file='top100genes_myel.csv')

#heatmaps
top5 <- mye.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(mye, features = top5$gene) + NoLegend()

DoHeatmap(mye, features = top5$gene) + scale_fill_gradientn(colors = something)

#####################################################################################


###########----------------COLORS--------------#################

#--> Previous Colors
#condition color = c("#1089eb", "#00ccb1", "#b52ec9", "#FF7F00", "#14690a", "#a62626")
NTF color = c("#1089eb", "#00ccb1", "#a62626")
#            c(alpha("#1089eb", 0.5), alpha("#00ccb1",0.5), alpha("#a62626",0.5))) 
#condition color = c("#1089eb", "#00ccb1", "#b52ec9", "#FF7F00", "#14690a", "#a62626")

magma = colorRampPalette(c('#050238','#59139e','#fa7916', '#fffb82', '#fffdbd'), alpha = TRUE)(8)
viridis = colorRampPalette(c('#fde725','#35b778','#3e4a89','#430154'), alpha = TRUE)(8)

magma = c("#050238FF", "#280963FF", "#4D108FFF", "#863077FF", "#CC5B3CFF", "#FA8B2DFF", "#FCC475FF", "#FFFDBD1A")
magma = c("#0502381A", "#280963FF", "#4D108FFF", "#863077FF", "#CC5B3CFF", "#FA8B2DFF", "#FCC475FF", "#FFFDBDFF")
magma = c("#0502381A", "#340B72FF", "#6F218AFF", "#CB5B3CFF", "#FB9E34FF", "#FEE872FF", "#FFFB9BFF", "#FFFDBDFF")


something = colorRampPalette(c('#000000','#000000','#33313b','#333333','#640E27','#900c3f','#ff8300','#ffd369', '#feffdb'), alpha = TRUE)(8)

NTF_colors = c('#f5bb06','#81007f','#e6194b')


# Load the "scales" package
require(scales)
# Create vector with levels of object@ident
identities <- unique(mye$PatientID)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(5)

################################################################


##########-----------Proportion barplot---------################

# --> estimate PDF size 
# --> save on export

# conditions --> c('seurat_clusters', 'NTF', 'PatientID')
library(dplyr)
library(tidyr)

meta <- mye@meta.data[,colnames(mye@meta.data) %in% c('seurat_clusters','PatientID')]
head(meta)

df <- meta %>% group_by(seurat_clusters,PatientID) %>% summarise(n = n()) %>% spread(PatientID, n, fill = 0)
df <- as.data.frame(df)
rownames(df) <- df$seurat_clusters
df <- df[, -which(names(df) %in% c("seurat_clusters"))]
df <- df/colSums(df)*100
df <- df/rowSums(df)*100

par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(as.matrix(t(df)), horiz=TRUE, col = my_color_palette, legend = colnames(df),  args.legend = list(x = "topright", bty = "n", inset = c(-0.25, 0)), border = FALSE)

head(mye@meta.data)

meta <- mye@meta.data[,colnames(mye@meta.data) %in% c('seurat_clusters','NTF')]
head(meta)

df <- meta %>% group_by(seurat_clusters,NTF) %>% summarise(n = n()) %>% spread(NTF, n, fill = 0)
df <- as.data.frame(df)
rownames(df) <- df$seurat_clusters
df <- df[, -which(names(df) %in% c("seurat_clusters"))]
df <- df/colSums(df)*100
df <- df/rowSums(df)*100

par(mfrow=c(1, 1), mar=c(3, 5, 1, 7))
barplot(as.matrix(t(df)), horiz=TRUE, col = NTF_colors , legend = colnames(df),  
        args.legend = list(x = "topright", bty = "n", inset = c(-0.23, 0.05)), border = FALSE)


################################################################

meta <- mye@meta.data[,colnames(mye@meta.data) %in% c('seurat_clusters','tech')]



##########----------------Saving plots--------------################

#setwd("/mnt/justinedata/HCC_fetal/figures/endo/")

# estimate pdf size on plots 
# change the width and height accordingly 

pdf(paste0("/mnt/justinedata/HCC_fetal/figures/mye/umapclsuters_11.pdf"), width=9, height=5.7)
DimPlot(mye, reduction = "umap", group.by = 'seurat_clusters',  pt.size=0.5)
dev.off()

genes = c('S100A8','FOLR2','SPP1','MT1G','CD163','DC1','DC2','pDC','C1QB')

print(paste0("saving expression plot --> /mnt/justinedata/HCC_fetal/figures/mye/"))
for (i in genes){
  print(i)
  pdf(paste0("/mnt/justinedata/HCC_fetal/figures/mye/expr_", toString(i),".pdf"), width=7, height=5)
  print(FeaturePlot(mye, min.cutoff = 0, pt.size = 0.5, features = i, 
                    blend.threshold = 0.01, cols = magma))
  dev.off()
  
}

print(paste0("saving violin plot --> /mnt/justinedata/HCC_fetal/figures/mye/"))
for (i in genes){
  print(i)
  pdf(paste0("/mnt/justinedata/HCC_fetal/figures/mye/vln_", toString(i),".pdf"), width=7, height=5)
  print(VlnPlot(mye, i, group.by = "seurat_clusters", pt.size = 0) + NoLegend())
  dev.off()
  
}


mye.markers <- FindAllMarkers(endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

print(paste0("saving heatmap --> /mnt/justinedata/HCC_fetal/figures/mye/"))
top5 <- mye.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pdf(paste0("/mnt/justinedata/HCC_fetal/figures/mye/heatmap.pdf"), width=11, height=8)
print(DoHeatmap(mye, features = top5$gene))
dev.off()

mye@meta.data["PatientID"] <- wholeatlas@meta.data[rownames(wholeatlas@meta.data) %in% rownames(mye@meta.data), "PatientID", drop = FALSE]
pdf(paste0("/mnt/justinedata/HCC_fetal/figures/mye/PatientID.pdf"), width=7, height=5)
DimPlot(mye, reduction = "umap", group.by = 'PatientID',  pt.size=0.5, cols = my_color_palette)
dev.off()

pdf(paste0("/mnt/justinedata/HCC_fetal/figures/mye/NTF.pdf"), width=9, height=5.7)
DimPlot(mye, reduction = "umap", group.by = 'NTF',cols = c('#81007f', '#f5bb06', '#e6194b'),  pt.size=0.5)
dev.off()

################################################################


##########----------------Saving plots--------------################

#setwd("/mnt/justinedata/HCC_fetal_norm_cirr/figures/endo/")

# estimate pdf size on plots 
# change the width and height accordingly 

pdf(paste0("/mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/umapcluster.pdf"), width=9, height=5.7)
DimPlot(mye, reduction = "umap", group.by = 'seurat_clusters',  pt.size=0.5)
dev.off()

pdf(paste0("/mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/umapcondition.pdf"), width=9, height=5.7)
DimPlot(mye, reduction = "umap", group.by = "condition", pt.size=0.5, cols = c('#f5bb06', "#f45905", '#81007f', "#9aceff",'#e6194b'))
dev.off()


genes = c('S100A8','FOLR2','SPP1','MT1G','CD163','DC1','DC2','pDC','C1QB')

print(paste0("saving expression plot --> /mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/"))
for (i in genes){
  print(i)
  pdf(paste0("/mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/expr_", toString(i),".pdf"), width=7, height=5)
  print(FeaturePlot(mye, min.cutoff = 0, pt.size = 0.5, features = i, 
                    blend.threshold = 0.01, cols = magma))
  dev.off()
  
}

print(paste0("saving violin plot --> /mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/"))
for (i in genes){
  print(i)
  pdf(paste0("/mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/vln_", toString(i),".pdf"), width=7, height=5)
  print(VlnPlot(mye, i, group.by = "seurat_clusters", pt.size = 0) + NoLegend())
  dev.off()
  
}


print(paste0("saving violin plot --> /mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/"))
for (i in genes){
  print(i)
  pdf(paste0("/mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/vln_", toString(i),".pdf"), width=7, height=5)
  print(VlnPlot(mye, i, group.by = "seurat_clusters", pt.size = 0) + NoLegend())
  dev.off()
  
}


genes = c('FOLR2','SPP1')
print(paste0("saving violin plot --> /mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/"))
for (i in genes){
  print(i)
  pdf(paste0("/mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/vln_", toString(i),"_TAMScond.pdf"), width=7, height=5)
  print(VlnPlot(mye[,mye$seurat_clusters %in% c(3,4,5,8,9) ], i, group.by = "condition", pt.size = 0) + NoLegend())
  dev.off()
  
}

mye.markers <- FindAllMarkers(mye, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

print(paste0("saving heatmap --> /mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/"))
top5 <- mye.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pdf(paste0("/mnt/justinedata/HCC_fetal_norm_cirr/figures/mye/heatmap.pdf"), width=11, height=10)
print(DoHeatmap(mye, features = top5$gene))
dev.off()



##########-----------Proportion barplot---------################

# --> estimate PDF size 
# --> save on export

# conditions --> c('seurat_clusters', 'NTF', 'PatientID')

meta <- mye@meta.data[,colnames(mye@meta.data) %in% c('seurat_clusters','condition')]
head(meta)

df <- meta %>% group_by(seurat_clusters,condition) %>% summarise(n = n()) %>% spread(condition, n, fill = 0)
df <- as.data.frame(df)
rownames(df) <- df$seurat_clusters
df <- df[, -which(names(df) %in% c("seurat_clusters"))]
df <- df/colSums(df)*100
df <- df/rowSums(df)*100

par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(as.matrix(t(df)), horiz=TRUE, col = my_color_palette, legend = colnames(df), cols = my_color_palette,
        args.legend = list(x = "topright", bty = "n", inset = c(-0.25, 0)), border = FALSE)

################################################################

VlnPlot(mye, c('CD209'), group.by = "condition", pt.size = 0.5) + NoLegend()

FeaturePlot(mye, min.cutoff = 0, pt.size = 0.5, features = 'CD209', 
            blend.threshold = 0.01, cols = magma)

table(mye[,mye$seurat_clusters %in% c(5,4,3)]$condition)



