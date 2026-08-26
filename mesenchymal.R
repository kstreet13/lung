################
# LOAD / RESET Mesenchymal SCE
################
require(HDF5Array)
require(SingleCellExperiment)
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
anno <- read.csv('data/umap_mes_coordinates.csv', row.names = 1)
anno$mes.clus <- paste0('M',anno$mes.clus)
sce <- sce[, match(rownames(anno), colnames(sce))]
all(colnames(sce) == rownames(anno)) # check
sce$clus.mes <- anno$mes.clus
reducedDim(sce,'umap.mes') <- anno[,c('UMAP1','UMAP2')]
rm(anno)
################


# UMAP plots
ind <- sample(ncol(sce))
plot(reducedDim(sce,'umap.mes')[ind,],asp=1,col=colorby(as.character(sce$clus.mes[ind])), cex=.5)
labelby(reducedDim(sce,'umap.mes'),as.character(sce$clus.mes))

plot(reducedDim(sce,'umap.mes')[ind,],asp=1,col=colorby(sce$condition)[ind], cex=.5)
legendby(sce$condition)

layout(matrix(1:2,nrow=1))
ind <- which(sce$condition=='RA')
plot(reducedDim(sce,'umap.mes')[ind,],asp=1,col=colorby(sce$condition)[ind], cex=.25)
ind <- which(sce$condition=='HO85')
plot(reducedDim(sce,'umap.mes')[ind,],asp=1,col=colorby(sce$condition)[ind], cex=.25)


# cluster x condition
barplot(table(sce$condition,sce$clus.mes), col=c(2,4), las=2, beside=TRUE)
barplot(table(sce$condition,sce$clus.mes), col=c(2,4), las=2, beside=FALSE)


# marker plots
gene <- 'Tgfb3'
plot(reducedDim(sce,'umap.mes')[ind,],asp=1,col=colorby(assay(sce,'binomial_deviance_residuals')[gene,ind], colors = c('grey90','grey90','lightgreen','green','darkgreen','blue','darkblue')), cex=.5, main=gene)

boxplot(assay(sce,'binomial_deviance_residuals')[gene,] ~ sce$clus.mes)

