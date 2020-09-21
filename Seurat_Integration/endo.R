
##########----------------READY? START!--------------################

endo <- readRDS('/mnt/justinedata/HCC_fetal/data/endo.RData')

endo <- subset(reference.integrated, idents = c('8','9','10','15','20','34'))
endo

DimPlot(endo, reduction = "umap", label = TRUE, pt.size = 2) + NoLegend()

##########----------------Usual Practice--------------################

endo <- RunPCA(endo, features = VariableFeatures(object = endo), npcs = 70)

print(endo[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(endo, dims = 1:2, reduction = "pca")

ElbowPlot(endo, ndims = 50)

endo <- FindNeighbors(endo, dims= 1:10)
endo <- FindClusters(endo, resolution = 0.6)
endo <- RunUMAP(endo, dims = 1:10, min.dist=0.001, spread=0.3)


saveRDS(endo, file = '/mnt/justinedata/HCC_fetal/data/endo.RData')

#####################################################################################



#############---------------------PLOTS-------------------################

cond = c('seurat_clusters', 'NTF', 'PatientID')

DimPlot(endo, reduction = "umap", label = TRUE, pt.size = 2, label.size = 6)# + NoLegend()

colorcode <- c('Adj Normal'= alpha("#1089eb", 0.3), 'fetal'= alpha("#00ccb1", 0.3), 'Tumor'= alpha("#a62626", 0.3))
DimPlot(endo, reduction = "umap", group.by = "NTF", pt.size=0.5,
        cols = colorcode)

genes = c('PLVAP','SPRY1','IGFBP3','TFF3','CLEC1B','SLC9A3R2','ACKR1')
FeaturePlot(endo, features = genes, min.cutoff = 0, pt.size = 0.1, blend.threshold = 0.01, 
            cols = magma, ncol = 4) 

genes = c('PLVAP','SPRY1','IGFBP3','TFF3','CLEC1B','SLC9A3R2','ACKR1','IL6')
FeaturePlot(endo, features = genes, min.cutoff = 0, pt.size = 0.1, blend.threshold = 0.01, 
            cols = magma, ncol = 4) 

VlnPlot(endo, genes, group.by = "seurat_clusters", pt.size = 0)

#####################################################################################



#############---------------------Renaming columns-------------------################

endometa <- endo@meta.data[,colnames(endo@meta.data) %in% c('seurat_clusters','NTF')]
head(endometa)

df <- endometa %>% 
  group_by(seurat_clusters,NTF) %>% 
  summarise(n = n()) %>% 
  spread(NTF, n, fill = 0)

df <- as.data.frame(df)

rownames(df) <- df$seurat_clusters

df <- df[,c('Adj Normal', 'fetal', 'Tumor')]

df <- df/colSums(df)*100
df <- df/rowSums(df)*100

par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(as.matrix(t(df)), horiz=TRUE, col=c("#1089eb", "#00ccb1", "#a62626"),
        legend = colnames(df),  
        args.legend = list(x = "topright", bty = "n", inset=c(-0.35, 0)))

#####################################################################################


##########-----------Saving plots---------################

#setwd("/mnt/justinedata/HCC_fetal/figures/endo/")

genes = c('PLVAP','SPRY1','IGFBP3','TFF3','CLEC1B','SLC9A3R2','ACKR1','HLA-DRA')

print(paste0("saving expression plot --> /mnt/justinedata/HCC_fetal/figures/endo/"))
for (i in genes){
  print(i)
  pdf(paste0("/mnt/justinedata/HCC_fetal/figures/endo/expr_", toString(i),".pdf"), width=7, height=5)
  print(FeaturePlot(endo, min.cutoff = 0, pt.size = 0.5, features = i, 
                    blend.threshold = 0.01, cols = magma))
  dev.off()
  
}

print(paste0("saving violin plot --> /mnt/justinedata/HCC_fetal/figures/endo/"))
for (i in genes){
  print(i)
  pdf(paste0("/mnt/justinedata/HCC_fetal/figures/endo/vln_", toString(i),".pdf"), width=7, height=5)
  print(VlnPlot(endo, i, group.by = "seurat_clusters", pt.size = 0) + NoLegend())
  dev.off()
  
}

endo@meta.data["PatientID"] <- wholeatlas@meta.data[rownames(wholeatlas@meta.data) %in% rownames(endo@meta.data), "PatientID", drop = FALSE]
pdf(paste0("/mnt/justinedata/HCC_fetal/figures/endo/PatientID.pdf"), width=9, height=5.7)
DimPlot(endo, reduction = "umap", group.by = 'PatientID',  pt.size=0.5)
dev.off()

pdf(paste0("/mnt/justinedata/HCC_fetal/figures/endo/NTF.pdf"), width=9, height=5.7)
DimPlot(endo, reduction = "umap", group.by = 'NTF',cols = c('#81007f', '#f5bb06', '#e6194b'),  pt.size=0.5)
dev.off()


endo.markers <- FindAllMarkers(endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

print(paste0("saving heatmap --> /mnt/justinedata/HCC_fetal/figures/endo/"))
top5 <- endo.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pdf(paste0("/mnt/justinedata/HCC_fetal/figures/endo/heatmap.pdf"), width=11, height=8)
print(DoHeatmap(endo, features = top5$gene))
dev.off()

################################################################


##########-----------Proportion barplot---------################

# --> estimate PDF size 
# --> save on export

# conditions --> c('seurat_clusters', 'NTF', 'PatientID')

meta <- endo@meta.data[,colnames(endo@meta.data) %in% c('seurat_clusters','PatientID')]
head(meta)

df <- meta %>% group_by(seurat_clusters,PatientID) %>% summarise(n = n()) %>% spread(PatientID, n, fill = 0)
df <- as.data.frame(df)
rownames(df) <- df$seurat_clusters
df <- df[, -which(names(df) %in% c("seurat_clusters"))]
df <- df/colSums(df)*100
df <- df/rowSums(df)*100

par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
barplot(as.matrix(t(df)), horiz=TRUE, col = my_color_palette, legend = colnames(df),  args.legend = list(x = "topright", bty = "n", inset = c(-0.25, 0)), border = FALSE)


meta <- endo@meta.data[,colnames(endo@meta.data) %in% c('seurat_clusters','NTF')]
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


###########----------------COLORS--------------#################

#--> Previous Colors
#condition color = c("#1089eb", "#00ccb1", "#b52ec9", "#FF7F00", "#14690a", "#a62626")
#NTF color = c("#1089eb", "#00ccb1", "#a62626")
#            c(alpha("#1089eb", 0.5), alpha("#00ccb1",0.5), alpha("#a62626",0.5))) 
#condition color = c("#1089eb", "#00ccb1", "#b52ec9", "#FF7F00", "#14690a", "#a62626")

magma = colorRampPalette(c('#050238','#59139e','#fa7916', '#fffdbd'), alpha = TRUE)(8)
viridis = colorRampPalette(c('#fde725','#35b778','#3e4a89','#430154'), alpha = TRUE)(8)

magma = c("#050238FF", "#280963FF", "#4D108FFF", "#863077FF", "#CC5B3CFF", "#FA8B2DFF", "#FCC475FF", "#FFFDBD1A")
magma = c("#0502381A", "#280963FF", "#4D108FFF", "#863077FF", "#CC5B3CFF", "#FA8B2DFF", "#FCC475FF", "#FFFDBDFF")

something = colorRampPalette(c('#000000','#000000','#33313b','#333333','#640E27','#900c3f','#ff8300','#ffd369', '#feffdb'), alpha = TRUE)(8)

NTF_color = c('#e6194b', '#f5bb06', '#81007f')

p <- DimPlot(endo, reduction = "umap", group.by = "PatientID", pt.size=0.5, label = TRUE)
patientid_color = unique(ggplot_build(p)$data[[1]]$colour)

################################################################



