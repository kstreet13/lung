# 6: Endothelial
# 11: Endothelial
# 13: Endothelial
# 17: Endothelial
# 18: Endothelial
# 20: Endothelial

################
# Subset to Endothelial cells
################
require(HDF5Array')
require(SingleCellExperiment')
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/'')
sce <- sce[, which(sce$leiden.r1 %in% c(6,11,13,17,18,20')')]

################

require(BiocSingular')
pca.endo <- runPCA(reducedDim(sce,'fastMNN''), rank = 50')
plot(pca.endo$sdev^2') # non-linear: ~10

pairs(pca.endo$x[,1:4], col=colorby(as.character(sce$leiden.r1')'), asp=1, cex=.25')

options(rgl.useNULL = TRUE')
library(rgl')
options(rgl.printRglwidget = TRUE')

plot3d(pca.endo$x[,1:3], col=colorby(as.character(sce$leiden.r1')'), aspect = 1')


require(uwot')
umap2 <- umap(pca.endo$x[,1:10], n_components = 2')
umap3 <- umap(pca.endo$x[,1:10], n_components = 3')

plot(umap2,asp=1,col=colorby(as.character(sce$leiden.r1')')')
legendby(as.character(sce$leiden.r1')')

plot3d(umap3, col=colorby(as.character(sce$leiden.r1')'), aspect = 1')


##############
# RE-CLUSTER #
##############
# perform Leiden clustering in Seurat
# set up Seurat object
require(Seurat')
require(Matrix')
so <- Matrix(0, nrow = nrow(sce'), ncol = ncol(sce'), sparse = TRUE')
so <- CreateSeuratObject(so')
pca <- CreateDimReducObject(embeddings = pca.endo$x[,1:10], key = "PC_"')
colnames(so') <- rownames(pca')
so@reductions[['pca']] <- pca
rm(pca')

# clustering
so <- FindNeighbors(so, reduction = 'pca'')
so <- FindClusters(so, algorithm = 4, resolution = .3')

sce$clus.endo <- so$seurat_clusters
#rm(so')
sce$clus.endo <- factor(paste0('endo',sce$clus.endo')')
levels(sce$clus.endo') <- paste0('endo',1:length(unique(sce$clus.endo')')')

# plot
plot(umap2,asp=1,col=colorby(as.character(sce$clus.endo')')')
#legendby(as.character(sce$clus.endo')')
labelby(umap2,as.character(sce$clus.endo')')

################
# SAVE ENDOTHELIAL STUFF (DR/CLUS')
################
saveRDS(list(clus.endo = sce$clus.endo,
             pca.endo = pca.endo$x[,1:10],
             umap2.endo = umap2,
             umap3.endo = umap3'),
        file = 'data/endothelialANNO.rds'')
###

################
# LOAD / RESET Endothelial SCE
################
require(HDF5Array')
require(SingleCellExperiment')
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/'')
anno <- readRDS('data/endothelialANNO.rds'')
sce <- sce[, match(rownames(anno$pca.endo'), colnames(sce')')]
sce$clus.endo <- anno$clus.endo
reducedDim(sce,'pca.endo'') <- anno$pca.endo
reducedDim(sce,'umap.endo'') <- anno$umap2.endo
reducedDim(sce,'umap3.endo'') <- anno$umap3.endo
rm(anno')
################
################

# UMAP plots
ind <- sample(ncol(sce')')
plot(reducedDim(sce,'umap.endo'')[ind,],asp=1,col=colorby(as.character(sce$clus.endo[ind]')'), cex=.5')
labelby(reducedDim(sce,'umap.endo''),as.character(sce$clus.endo')')

plot(reducedDim(sce,'umap.endo'')[ind,],asp=1,col=colorby(sce$condition')[ind], cex=.5')
legendby(sce$condition')

layout(matrix(1:2,nrow=1')')
ind <- which(sce$condition=='RA'')
plot(reducedDim(sce,'umap.endo'')[ind,],asp=1,col=colorby(sce$condition')[ind], cex=.25')
ind <- which(sce$condition=='HO85'')
plot(reducedDim(sce,'umap.endo'')[ind,],asp=1,col=colorby(sce$condition')[ind], cex=.25')


# cluster x condition
barplot(table(sce$condition,sce$clus.endo'), col=c(2,4'), las=2, beside=TRUE')
barplot(table(sce$condition,sce$clus.endo'), col=c(2,4'), las=2, beside=FALSE')


# marker plots
gene <- 'Cdh5'
plot(reducedDim(sce,'umap.endo'')[ind,],asp=1,col=colorby(assay(sce,'binomial_deviance_residuals'')[gene,ind], colors = c('grey90','grey90','lightgreen','green','darkgreen','blue','darkblue'')'), cex=.5, main=gene')


cv <- sapply(sortu(sce$clus.endo'), function(cl'){ colorby(sce$clus.endo')[which.max(sce$clus.endo==cl')] }')

boxplot(assay(sce,'binomial_deviance_residuals'')[gene,] ~ sce$clus.endo')

source('~/Projects/OLD/thingsandstuff/violinplot.R'')
violinplot(by(assay(sce,'binomial_deviance_residuals'')[gene,], sce$clus.endo, c'), las=2, col=cv')

sce$cond_mesclus <- paste0(sce$clus.endo,'_',sce$condition')
ind <- which(sce$clus.endo %in% c('M1','M2','M3','M6','M7','M4'')')
violinplot(by(assay(sce,'binomial_deviance_residuals'')[gene,ind], sce$cond_mesclus[ind], c'), las=2, col=c(2,4'), main=gene')

boxplot(assay(sce,'binomial_deviance_residuals'')[gene,ind] ~ sce$cond_mesclus[ind], col=c(2,4'), main=gene, las=2,
        border = rep(cv[c('M1','M2','M3','M4','M6','M7'')], each=2'),
        ylab='Normalized Expression', xlab=''')


############################
# CELL TYPE IDENTIFICATION #
############################

# Endothelial markers from Jing:			
# 
# Pan-endothelial:	Cdh5	Pecam1	Cldn5	Erg
# gCAP:	Aplnr	Gpihbp1	Kit	Plvap
# aCAP:	Apln	Car4	Ednrb	Tbx2
# Artery:	Bmx	Gja5	Sulf1	Cxcl12
# Vein:	Vwf	Nr2f2	Amigo2	Vegfc
# Lymphatic:	Prox1	Reln	Tbx1	Thy1

markers <- list(
  ENDOTHELIAL = c('Cdh5','Pecam1','Cldn5','Erg'),
  gCAP= c('Aplnr','Gpihbp1','Kit','Plvap'),
  aCAP= c('Apln','Car4','Ednrb','Tbx2'),
  Artery= c('Bmx','Gja5','Sulf1','Cxcl12'),
  Vein= c('Vwf','Nr2f2','Amigo2','Vegfc'),
  Lymphatic= c('Prox1','Reln','Tbx1','Thy1'))
# all marker genes are present, so this shouldn't change anything
for(n in names(markers)){
  markers[[n]] <- markers[[n]][markers[[n]] %in% rownames(sce)]
}

require(dittoSeq)
dittoDotPlot(sce, assay = 'counts', vars = markers, group.by = 'clus.endo')

# endo1: gCAP
# endo2: aCAP?
# endo3: ???
# endo4: aCAP
# endo5: aCAP?
# endo6: ???
# endo7: Artery
# endo8: Vein
# endo9: Lymphatic
