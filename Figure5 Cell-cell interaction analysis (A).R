#-------------------------NICHENET-----------------------###
load('hepa-endo_NicheNet.RData')

library(nichenetr)
library(Seurat)
library(tidyverse)
library(ggplot2)

nichenet_output$ligand_activity_target_heatmap

###------------------------------------------------------###

###----------------------CirclePlot----------------------###
library(circlize)

pval<-read.csv('./data/pvalues.txt', sep='\t', row.names=1, check.names=FALSE, header = TRUE)
pval <- pval[,!(colnames(pval) %in% colnames(pval)[0:8])]

lsname <- c()
for (i in colnames(pval)){
  if (strsplit(i, '_')[[1]][2] == strsplit(i, '_')[[1]][4]){
    lsname <- c(lsname,i)
  }
}

pval <- pval[,lsname]
pval[pval == 0] <- 0.001
pval <- abs(log10(pval))
max(pval);min(pval)
pval[pval < 2] <- 0
apply(pval, 2, function(x)table(x))

df <- apply(pval, 2, function(x)sum(x>=2))
table(df)
df = df[df > 0]
table(df)

Ndf <- data.frame(from=c(), to=c(), Freq=c())
Tdf <- data.frame(from=c(), to=c(), Freq=c())

for (i in names(df)){
  if (strsplit(i, '_')[[1]][2] == strsplit(i, '_')[[1]][4]){
    from = strsplit(i, '_')[[1]][1]
    to   = strsplit(i, '_')[[1]][3]
    if(strsplit(i, '_')[[1]][2] == 'Normal'){
      Ndf <- rbind(Ndf,data.frame(from = from, to = to, Freq=df[[i]]))
    }else if(strsplit(i, '_')[[1]][2] == 'Tumor'){
      Tdf <- rbind(Tdf,data.frame(from = from, to = to, Freq=df[[i]]))
    }else{
      print('Weird Stuff')
    }
    
  }
}

Ncol = c(CD4   = '#293462',
         CD8   = '#a64942',
         NK    = '#8a00d4',
         dc    = '#d527b7',
         endo  = '#fe5f55',
         fibro = '#9ea9f0',
         hepa  = '#ffc6be',
         mac   = '#fff1c1')

Tcol = c(CD4   = '#0c084c',
         CD8   = '#096386',
         NK    = '#616f39',
         dc    = '#76a21e',
         endo  = '#6fc2d0',
         fibro = '#beeef7',
         hepa  = '#f6e0b3',
         mac   = '#fab95b')

chordDiagram(Ndf, grid.col = Ncol, title("Normal", cex = 0.8), transparency = 0.3)
chordDiagram(Tdf, grid.col = Tcol, title("Tumor", cex = 0.8), transparency = 0.3)

###------------------------------------------------------###
