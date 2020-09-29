library(corrplot)
library(stringr)
library(RColorBrewer)

CMM <- read.csv('CMM_binaryMat.csv', check.names=FALSE, row.names=1)
class(CMM)
CMM[1:5,1:5]
rownames(CMM)[1:5] 
ll <- sapply(rownames(CMM) ,function(x) {str_split(x , '[_ ]')[[1]][1]})
rn = as.data.frame(rownames(CMM),ll)
rn['short']<- ll

df <- read.csv('./data/HMD_HumanPhenotype.rpt.txt', sep = '\t', header = FALSE)
df <- df[,c(1,5)]
colnames(df) <- c('human','mouse')
head(df)
dfgenes <- df[df$mouse %in% ll,]
df1 <- merge(rn, dfgenes, by.x = 'short', by.y = 'mouse', all = TRUE)

CMM <- CMM[!is.na(df1$human),]
df1 <- df1[!is.na(df1$human),]
rownames(CMM) <- df1$human
CMM['mergerow']<-rownames(CMM)
CMM$mergerow

HCCFM <- read.csv('HCCF_binaryMat.csv', check.names=FALSE, row.names=1)
class(HCCFM)
HCCFM[1:5,1:5]
rownames(HCCFM)[1:5] 
ll <- sapply(rownames(HCCFM) ,function(x) {str_split(x , '[_ ]')[[1]][1]})
rn = as.data.frame(rownames(HCCFM),ll)
rn['short']<- ll
rownames(HCCFM) <- rn$short
HCCFM['mergerow']<-rownames(HCCFM)
HCCFM$mergerow

mat <- merge(CMM, HCCFM, by='mergerow', all = TRUE)
rownames(mat) <- mat$mergerow
mat <- mat[, !colnames(mat) %in% c("mergerow")]
mat[is.na(mat)] <- 0
mat[1:5,1:5]

CMM <- CMM[, !colnames(CMM) %in% c("mergerow")]
HCCFM <- HCCFM[, !colnames(HCCFM) %in% c("mergerow")]


##----Correlation between Transcription Factors----##

regulon <- t(mat)
regulon[1:5,1:5]
class(regulon)

corr_mat_reg=cor(regulon,method="s")
corr_mat_reg[1:5,1:5]
ord=hclust(1-as.dist(corr_mat_reg))$order

if(TRUE){
  pdf('corrplt3.pdf', height = 30, width = 30)
  corrplot(corr_mat_reg[ord,ord], title = "Correlation Plot", method = "color", outline = T, 
           addgrid.col = "white", order="hclust", mar = c(4,0,4,0), addrect = 8, 
           rect.col = "black", rect.lwd = 5,cl.pos = "b",
           col = colorRampPalette(c('#ffffff','#ffffff','#fcfcd2', '#ffffc4','#ffffb8',
                                    '#ffffa3','#ffffa6', "#C7E9B4", "#7FCDBB","#41B6C4", 
                                    "#1D91C0", "#225EA8" ,"#253494" ,"#081D58"))(100))
  
  dev.off()
}



##----Correlation of orthologues gene expression----##

df <- as.data.frame(x=colnames(CMM), colnames(CMM))
colnames(df) <- 'ID' 
df['categories']<- 'mouse'

df1 <- as.data.frame(colnames(HCCFM), colnames(HCCFM))
colnames(df1) <- 'ID' 
df1['categories']<- 'human'
library(plyr)

info <- rbind.fill(df, df1)

mat1 <- as.data.frame(t(mat))
mat1['rn']<-rownames(mat1)

df <- merge(mat1, info, by.x = 'rn', by.y = 'ID', all=TRUE)
rownames(df) <- df$rn
head(df)
mat1 <- df[, !colnames(df) %in% c("rn")]
mat1['categories'] <- as.factor(mat1$categories)
table(mat1$categories)
mat1 <- mat1[, !colnames(mat1) %in% c("categories")]

mat1 <- mat1[rownames(mat1) %in% rownames(cellInfo),]
downsam <- t(mat1[sample(nrow(mat1), 800), ])
downsam <- t(mat1[rownames(mat1) %in% rownames(cellInfo), ])

downsam = cor(downsam,method="pearson")
downsam[is.na(downsam)] <- 0
downsam[1:5,1:5]
dim(downsam)
ord=hclust(1-as.dist(downsam))$order
downsam <- downsam[ord,ord]

pdf('CTcorrplt.pdf', height = 50, width = 50)
corrplot(downsam[ord,ord], title = "Correlation Plot", method = "color", outline = T, 
         addgrid.col = "white", order="hclust", mar = c(4,0,4,0), addrect = 8, 
         rect.col = "black", rect.lwd = 5,cl.pos = "b",
         col = colorRampPalette(c('#ffffff','#ffffff','#fcfcd2', '#ffffc4','#ffffb8',
                                  '#ffffa3','#ffffa6', "#C7E9B4", "#7FCDBB","#41B6C4", 
                                  "#1D91C0", "#225EA8" ,"#253494" ,"#081D58"))(100))

dev.off()

##--------------------------------------------##

library(reshape2)
melted_cormat <- melt(downsam)
head(melted_cormat)

library(ggplot2)
library(viridis)
library(pheatmap)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_viridis(discrete=FALSE)


info1 <- info[info$ID %in% rownames(downsam),'categories', drop=FALSE]

cellInfo <- read.csv('CellInfo.csv', check.names=FALSE, row.names=1)
head(cellInfo)
cellInfo['ID'] <- rownames(cellInfo)
cellInfo <- merge(cellInfo, info, by = 'ID', all=TRUE)
cellInfo <- cellInfo[cellInfo$CellType %in% c('Tom+ Mac','Tom- Mac','FLM','TAM1','TAM3','TAM2'),]
unique(cellInfo$CellType)
rownames(cellInfo)<-cellInfo$ID
cellInfo <- cellInfo[,!colnames(cellInfo) %in% 'tree']
colnames(cellInfo) <- c('ID','CellType','Species')
cellInfo <- cellInfo[rownames(cellInfo) %in% rownames(df),]

levels(cellInfo$CellType)<-c(levels(cellInfo$CellType),'Tompos_Mac', "Tomneg_Mac") 
cellInfo[cellInfo == "Tom+ Mac"] <- NA
cellInfo[is.na(cellInfo)] <- 'Tompos_Mac'
cellInfo[cellInfo == "Tom- Mac"] <- NA
cellInfo[is.na(cellInfo)] <- "Tomneg_Mac"

unique(cellInfo$Species)
unique(cellInfo$CellType)

infoCT <- cellInfo

my_colour = list(
   Species = c(human = '#1f77b4', mouse = '#ff7f0e'),
   CellType = c(Tompos_Mac = '#e377c2', Tomneg_Mac = '#279e68', FLM = '#7d010f',
                 TAM1 = '#d62728', TAM3 = '#8c564b', TAM2 = '#aa40fc')
)

pdf('CTcorrplttest.pdf', height = 50, width = 50)
pheatmap(downsam[ord,ord], annotation_colors = my_colour, annotation_col = cellInfo, 
         annotation_row = cellInfo, show_rownames=F, show_colnames=F)
dev.off()


p <- pheatmap(downsam, annotation_colors = my_colour, annotation_col = cellInfo, 
              annotation_row = cellInfo, show_rownames=F, show_colnames=F)

plot(p$tree_row)
abline(h=11.2, col="red", lty=2, lwd=2)
table(p$tree_row$labels %in% rownames(cellInfo) )

df <- as.data.frame(cutree(p$tree_row, h=11.2))
colnames(df) <-'tree'
df['ID'] <- rownames(df)
table(df$ID %in% rownames(cellInfo))
cellInfo['ID'] <- rownames(cellInfo)
cellInfo <- merge(cellInfo, df, by = 'ID', all=TRUE)
rownames(cellInfo) <- cellInfo$ID

cellInfo <- cellInfo[,!colnames(cellInfo) %in% 'ID']
cellInfo['tree'] = as.factor(cellInfo$tree)

pdf('CTcorrplttest.pdf', height = 50, width = 50)
p <- pheatmap(downsam, annotation_colors = my_colour, annotation_col = cellInfo, 
              annotation_row = cellInfo, show_rownames=F, show_colnames=F)
dev.off()
