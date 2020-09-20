setwd('/Users/justineseowjiawen/Dropbox/Liver_Paper/nanostring/')

library(pheatmap)
library(ggbiplot)
library(reshape2)

SNR <- read.csv('data/comb_SNR_wFetal.csv', header = TRUE, row.names=1)
cellinfo <- read.csv('data/sampleinfo_wFetal.csv', header = TRUE, row.names=1)

df <- SNR
df <- df[,!(colnames(df) %in% c('P14_N_S7'))]
cellinfo <- cellinfo[!(rownames(cellinfo) %in% c('P14_N_S7')),]

if(TRUE) {
  pca <- prcomp(t(dist(df)))
  rot <- data.frame(pca[["rotation"]])
  coor <- rot[,colnames(rot) %in% c('PC2','PC3')]
  pheatmap(coor, display_numbers = TRUE, main = 'SNR Normalization')
} 
if(TRUE) {
  plot(coor$PC2, coor$PC3, cex=0)
  text(coor$PC2, coor$PC3, labels=rownames(coor), cex = 0.5)
}

SNRhighvar <- c('AFP', 'UBB', 'CD74', 'H3F3A', 'HLA−E', 'HLA−DRB', 'ARG1', 'OAZ1', 'HES1', 'TIGIT', 
                'ALB', 'LY6E', 'pan−Melanocyte', 'STAT1', 'STAT3', 'PTEN', 'HIF1A', 'ICAM1', 'BATF3', 
                'CXCL10', 'STAT2', 'FOLR2', 'VEGFR2', 'CTNNB1', 'DLL4', 'PLVAP', 'AKT1', 'SDHA',
                'CD68', 'CD4', 'Notch2', 'VEGFA', 'RAB7A', 'SPP1', 'IFNG', 'CSF1R', 'GZMB', 'LAG3', 
                'NKG7', 'PDL1', 'VISTA', 'EPCAM','CD3E', 'CMKLR1', 'PDL2', 'CXCR6', 'CTLA4', 'IL12B', 
                'Tim3', 'CD11b', 'DKK2', 'FAS', 'FOXP3', 'IL15', 'IL6')

df1 <- df[c(rownames(df) %in% SNRhighvar),]

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

ann_colors = list(
  NormalvsTumor = c(Normal = '#5000ff', Tumor = '#e6194B', Fetal = '#d6b8f5'),
  Patient_ID = c(P8 = '#17becf', P14 = '#c5b0d5', P15 = '#c49c94', F14 = '#ffbb78', F14_5 = '#98df8a', F16_4 = '#ff9896')
)
col <- colorRampPalette(c('#000000',"#1818c4","#39c439","#fff129","#d10f0f",
                          "#8f39c4","#40025c"))(60)

##------------Heatmap------------##
p <- pheatmap(as.matrix(df1), annotation_col = cellinfo, angle_col = "45", 
              scale = "row", clustering_distance_rows = "manhattan", 
              annotation_colors = ann_colors, clustering_callback = callback,
              legend = TRUE, fontsize_row = 6, fontsize_col = 6, color = col, border_color = FALSE, 
              main = 'SNR heatmap')

##------------Top barplot------------##
df2 <- as.data.frame(t(apply(df1,1, scale)))
colnames(df2) <- colnames(df1)
CT <- as.data.frame(t(df2[c('FOLR2','PLVAP','VEGFA', 'VEGFR2', 'Notch2', 'DLL4', 'FOXP3', 'CTLA4'),p$tree_col$order]))
CT['sampleid'] <- rownames(CT)
CT2 <- melt(CT)
colnames(CT2) <- c('SampleID','Gene','Value')
ggplot(CT2, aes(fill=Gene, y=Value, x=SampleID)) + 
  geom_bar(position="stack", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(limits = rownames(CT))  


##------------Proportional barplot------------##
cellinfo <- cellinfo[,colnames(cellinfo) %in% c('Patient_ID','NormalvsTumor')]

plot(p$tree_col)
abline(h=25, col="red", lty=2, lwd=2)
table(p$tree_col$labels %in% rownames(cellinfo) )

df <- as.data.frame(cutree(p$tree_col, h=25))
colnames(df) <-'tree'
df['ID'] <- rownames(df)
table(df$ID %in% rownames(cellinfo))
cellinfo['ID'] <- rownames(cellinfo)
cellinfo <- merge(cellinfo, df, by = 'ID', all=TRUE)
rownames(cellinfo) <- cellinfo$ID

cellinfo$tree <- as.character(cellinfo$tree)
cellinfo$tree[cellinfo$tree == "1"] <- "OncoFetal Ecosystem"
cellinfo$tree[cellinfo$tree == "2"] <- "Immune Exclusion"
cellinfo['condition'] = as.factor(cellinfo$tree)

cellinfo$NormalvsTumor <- as.character(cellinfo$NormalvsTumor)
cellinfo$NormalvsTumor[cellinfo$NormalvsTumor == "Normal"]  <- 'Adj. Normal'
cellinfo['NormalvsTumor'] = as.factor(cellinfo$NormalvsTumor)
cellinfo$NormalvsTumor <- factor(cellinfo$NormalvsTumor, 
                    levels=c("Fetal","Tumor","Adj. Normal" ))
cellinfo$condition <- factor(cellinfo$condition, 
                                 levels=c("OncoFetal Ecosystem", "Immune Exclusion"))

library(dplyr)
library(tidyr)

meta <- cellinfo[,colnames(cellinfo) %in% c('condition','NormalvsTumor')]
head(meta)

df <- meta %>% group_by(condition,NormalvsTumor) %>% summarise(n = n()) %>% spread(NormalvsTumor, n, fill = 0)
df <- as.data.frame(df)
rownames(df) <- df$condition
df <- df[, -which(names(df) %in% c("condition"))]
df <- df/colSums(df)*100
df <- df/rowSums(df)*100

NormalvsTumor = c('#d6b8f5','#e6194B','#5000ff')
par(mfrow=c(1, 1), mar=c(5, 3, 4, 8))
barplot(as.matrix(t(df)), horiz=FALSE, legend = colnames(df),  col=NormalvsTumor, 
        args.legend = list(x = "topright", bty = "n", inset = c(-0.3, 0)), border = FALSE)

#####-----------------------------------------------------------------------####









