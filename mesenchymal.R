################
# LOAD / RESET Mesenchymal SCE
################
require(HDF5Array)
require(SingleCellExperiment)
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
anno <- read.csv('data/umap_mes_coordinates.csv', row.names = 1)
#anno$mes.clus <- paste0('M',anno$mes.clus)
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
gene <- 'Thbs1'
plot(reducedDim(sce,'umap.mes')[ind,],asp=1,col=colorby(assay(sce,'binomial_deviance_residuals')[gene,ind], colors = c('grey90','grey90','lightgreen','green','darkgreen','blue','darkblue')), cex=.5, main=gene)


cv <- sapply(sortu(sce$clus.mes), function(cl){ colorby(sce$clus.mes)[which.max(sce$clus.mes==cl)] })

boxplot(assay(sce,'binomial_deviance_residuals')[gene,] ~ sce$clus.mes)

source('~/Projects/OLD/thingsandstuff/violinplot.R')
violinplot(by(assay(sce,'binomial_deviance_residuals')[gene,], sce$clus.mes, c), las=2, col=cv)

sce$cond_mesclus <- paste0(sce$clus.mes,'_',sce$condition)
ind <- which(sce$clus.mes %in% c('M1','M2','M3','M6','M7','M4'))
violinplot(by(assay(sce,'binomial_deviance_residuals')[gene,ind], sce$cond_mesclus[ind], c), las=2, col=c(2,4), main=gene)

boxplot(assay(sce,'binomial_deviance_residuals')[gene,ind] ~ sce$cond_mesclus[ind], col=c(2,4), main=gene, las=2,
        border = rep(cv[c('M1','M2','M3','M4','M6','M7')], each=2),
        ylab='Normalized Expression', xlab='')


