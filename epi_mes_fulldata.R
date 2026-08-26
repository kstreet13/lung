################
# Epithelial/Mesenchymal clusters in full dataset context
################
require(HDF5Array)
require(SingleCellExperiment)

# Epithelial setup
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
sce <- sce[ ,which(sce$leiden.r1 %in% c(2,14,21,22,32,33,35))]
anno <- readRDS('data/epithelialANNO.rds')
sce$clus.epi <- anno$clus.epi
sce$epi.col <- colorby(sce$clus.epi)
rm(anno)
epi <- sce

# Mesenchymal setup
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
anno <- read.csv('data/umap_mes_coordinates.csv', row.names = 1)
anno$mes.clus <- paste0('M',anno$mes.clus)
sce <- sce[, match(rownames(anno), colnames(sce))]
all(colnames(sce) == rownames(anno)) # check
sce$clus.mes <- anno$mes.clus
reducedDim(sce,'umap.mes') <- anno[,c('UMAP1','UMAP2')]
rm(anno)
mes <- sce

# LOAD / RESET FULL SCE
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')

# INTEGRATE
# epi
df <- colData(epi)[,c('clus.epi','epi.col')]
colData(sce) <- cbind(colData(sce), df[match(colnames(sce),rownames(df)),])
# mes
df <- colData(mes)[,c('clus.mes'), drop=FALSE]
colData(sce) <- cbind(colData(sce), df[match(colnames(sce),rownames(df)), ,drop=FALSE])

sce$clus.epi <- as.character(sce$clus.epi)
rm(epi,mes,df)
################


################
# PLOTS
################
ind <- sample(ncol(sce))

# UMAP by sample
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = colorby(as.character(sce$sample[ind])), xlab = 'UMAP-1', ylab = 'UMAP-2')
legendby(as.character(sce$sample))

# UMAP by condition
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = colorby(as.character(sce$condition[ind])), xlab = 'UMAP-1', ylab = 'UMAP-2')
legendby(as.character(sce$condition))

# UMAP by cluster
cc <- rep(c(brewer.pal(9,'Set1'), brewer.pal(8,'Set2'), brewer.pal(12,'Set3')), length.out = length(levels(sce$leiden.r1)))
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = cc[sce$leiden.r1[ind]], xlab = 'UMAP-1', ylab = 'UMAP-2')

labelby_small <- function(coords, x){
  x <- factor(x)
  if(length(levels(x)) == length(x)){
    warning('labelby: x contains only unique values')
  }
  cc.full <- colorby(x)
  cc <- sapply(levels(x), function(clID){
    cc.full[which.max(x==clID)]
  })
  
  centers <- t(sapply(levels(x), function(clID){
    colMedians(coords[which(x==clID),])
  }))
  points(centers,pch=1,cex=1.5)
  points(centers,pch=16,cex=1.5, col=cc)
  text(centers, labels = levels(x), col = 1, font=2, cex=.75)
}

# Epi/Mes clusters
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = 'grey80', xlab = 'UMAP-1', ylab = 'UMAP-2')

points(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = colorby(sce$clus.epi[ind]), xlab = 'UMAP-1', ylab = 'UMAP-2')
labelby_small(reducedDim(sce,'umap'), sce$clus.epi)

points(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = colorby(sce$clus.mes[ind]), xlab = 'UMAP-1', ylab = 'UMAP-2')
labelby_small(reducedDim(sce,'umap'), sce$clus.mes)


# marker plots
gene <- 'Thbs1'
plot(reducedDim(sce,'umap')[ind,],asp=1,col=colorby(assay(sce,'binomial_deviance_residuals')[gene,ind], colors = c('grey90','lightgreen','green','darkgreen','blue','darkblue')), cex=.25)


