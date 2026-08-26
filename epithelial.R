# examine epithelial cells, especially AT1/AT2
# 2: AT2
# 14: AT2
# 21: AT1
# 22: AT2
# 32: AT2
# 33: AT1
# 35: Club/Ciliated

set.seed(1)

require(HDF5Array)
require(SingleCellExperiment)
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
sce <- sce[ ,which(sce$leiden.r1 %in% c(2,14,21,22,32,33,35))]

require(BiocSingular)
pca.epi <- runPCA(reducedDim(sce,'fastMNN'), rank = 50)
plot(pca.epi$sdev^2) # non-linear: ~12

pairs(pca.epi$x[,1:4], col=colorby(as.character(sce$leiden.r1)), asp=1, cex=.25)

options(rgl.useNULL = TRUE)
library(rgl)
options(rgl.printRglwidget = TRUE)

plot3d(pca.epi$x[,1:3], col=colorby(as.character(sce$leiden.r1)), aspect = 1)


require(uwot)
umap2 <- umap(pca.epi$x[,1:12], n_components = 2)
umap3 <- umap(pca.epi$x[,1:12], n_components = 3)

plot(umap2,asp=1,col=colorby(as.character(sce$leiden.r1)))
legendby(as.character(sce$leiden.r1))

plot3d(umap3, col=colorby(as.character(sce$leiden.r1)), aspect = 1)


# multiple umaps to examine common features
# umaps <- lapply(c(8,14,20,50), function(k){
#   umap <- umap(pca$x[,1:k], n_components = 2)
#   return(umap)
# })
# layout(matrix(1:4, 2,2, byrow = TRUE))
# for(umap in umaps){
#   ind <- sample(nrow(umap))
#   plot(umap[ind,],asp=1,col=colorby(as.character(sce$condition[ind])), cex=.25)
#   #legendby(as.character(sce$leiden.r1))
# }
# layout(1)

##############
# RE-CLUSTER #
##############
# perform Leiden clustering in Seurat
# set up Seurat object
require(Seurat)
require(Matrix)
so <- Matrix(0, nrow = nrow(sce), ncol = ncol(sce), sparse = TRUE)
so <- CreateSeuratObject(so)
pca <- CreateDimReducObject(embeddings = pca.epi$x[,1:12], key = "PC_")
colnames(so) <- rownames(pca)
so@reductions[['pca']] <- pca
rm(pca)

# clustering
so <- FindNeighbors(so, reduction = 'pca')
so <- FindClusters(so, algorithm = 4, resolution = .4)

sce$clus.epi <- so$seurat_clusters
#rm(so)
sce$clus.epi <- factor(paste0('E',sce$clus.epi))
levels(sce$clus.epi) <- paste0('E',1:11)

# plot
plot(umap2,asp=1,col=colorby(as.character(sce$clus.epi)))
#legendby(as.character(sce$clus.epi))
labelby(umap2,as.character(sce$clus.epi))

plot3d(umap3, col=colorby(as.character(sce$clus.epi)), aspect = 1)

################
# SAVE EPITHELIAL STUFF (DR/CLUS)
################
saveRDS(list(clus.epi = sce$clus.epi,
             pca.epi = pca.epi$x[,1:12],
             umap2.epi = umap2,
             umap3.epi = umap3),
        file = 'data/epithelialANNO.rds')
###


################
# LOAD / RESET Epithelial SCE
################
require(HDF5Array)
require(SingleCellExperiment)
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
sce <- sce[ ,which(sce$leiden.r1 %in% c(2,14,21,22,32,33,35))]
anno <- readRDS('data/epithelialANNO.rds')
sce$clus.epi <- anno$clus.epi
reducedDim(sce,'pca.epi') <- anno$pca.epi
reducedDim(sce,'umap2.epi') <- anno$umap2.epi
reducedDim(sce,'umap3.epi') <- anno$umap3.epi
rm(anno)
################
################

# UMAP plots
ind <- sample(ncol(sce))
plot(reducedDim(sce,'umap2.epi')[ind,],asp=1,col=colorby(as.character(sce$clus.epi[ind])), cex=.5)
labelby(reducedDim(sce,'umap2.epi'),as.character(sce$clus.epi))

plot(reducedDim(sce,'umap2.epi')[ind,],asp=1,col=colorby(sce$condition)[ind], cex=.5)
legendby(sce$condition)


layout(matrix(1:2,nrow=1))
ind <- which(sce$condition=='RA')
plot(reducedDim(sce,'umap2.epi')[ind,],asp=1,col=colorby(sce$condition)[ind], cex=.25)
ind <- which(sce$condition=='HO85')
plot(reducedDim(sce,'umap2.epi')[ind,],asp=1,col=colorby(sce$condition)[ind], cex=.25)



plot3d(reducedDim(sce,'pca.epi')[,1:3], col=colorby(as.character(sce$clus.epi)), aspect = 1)

# marker plots
gene <- 'Tgfbr2'
plot(reducedDim(sce,'umap2.epi')[ind,],asp=1,col=colorby(assay(sce,'binomial_deviance_residuals')[gene,ind], colors = c('grey90','lightgreen','green','darkgreen','blue','darkblue')), cex=.5, main=gene)

# "Tgfbr2"   "Tgfbr3"   "Tgfb2"    "Tgfbr1"   "Tgfbi"   "Tgfb1i1"  "Tgfb1"    "Tgfbrap1" "Tgfb3"    "Tgfbr3l" 

boxplot(assay(sce,'binomial_deviance_residuals')[gene,] ~ sce$clus.epi)


# cluster x condition barplot
barplot(table(sce$condition,sce$clus.mes), col=c(2,4), las=2, beside=TRUE)
barplot(table(sce$condition,sce$clus.mes), col=c(2,4), las=2, beside=FALSE)


# markers from LungMAP
# AT2 (mouse):
# Lamp3 Abca3 Kcnj15 Sftpa1 Lyz2 Lyz3 Lyz1 Sftpb Sftpc Slc34a2 Chi3l1
# AT1 (mouse)
# Ager Rtkn2 Sema3b Akap5 Cldn18 Emp2 Aqp5 Clic5 Msln Lmo7 Hopx
# proliferation markers
# Top2a, Mki67, Pcna

markers <- list(
  AT1 = c('Ager','Rtkn2','Sema3b','Akap5','Cldn18','Emp2','Aqp5','Clic5','Msln','Lmo7','Hopx'),
  AT2 = c('Lamp3','Abca3','Kcnj15','Sftpa1','Lyz2','Lyz3','Lyz1','Sftpb','Sftpc','Slc34a2','Chi3l1'),
  ClubCil = c('Scgb1a1','Scgb3a2','Foxj1','Dynlrb2'),
  Other = c('Retnla')
)
markers$AT1 <- markers$AT1[markers$AT1 %in% rownames(sce)]
markers$AT2 <- markers$AT2[markers$AT2 %in% rownames(sce)]
markers$ClubCil <- markers$ClubCil[markers$ClubCil %in% rownames(sce)]

sce$clus_cond <- factor(paste0(sce$clus.epi,'_',sce$condition))

require(dittoSeq)
dittoDotPlot(sce, assay = 'counts', vars = markers, group.by = 'clus.epi')

# E9 (pale orange) is club/ciliated


##################
# AT1/AT2 SCORES #
##################
so <- sce
assayNames(so)[2] <- 'logcounts'
assay(so,'counts') <- as(assay(so,'counts'), 'dgCMatrix')
assay(so,'logcounts') <- as.matrix(assay(so,'logcounts'))
so <- as.Seurat(so)

so <- AddModuleScore(so, features = markers, name = c('AT1_','AT2_','ClubCil'))

layout(matrix(1:2,nrow=1))
plot(umap2,asp=1,col=colorby(so$AT1_1, colors = c('grey80','yellow', 'green','blue','darkblue')), main = 'AT1 Score', cex=.5)
plot(umap2,asp=1,col=colorby(so$AT2_2, colors = c('grey80','yellow', 'green','blue','darkblue')), main = 'AT2 Score', cex=.5)
layout(1)

# difference
plot(umap2,asp=1,col=colorby(so$AT2_2-so$AT1_1, colors = c('darkred','red', 'grey90','blue','darkblue')), main = 'AT2 - AT1 Score')
ind <- sample(ncol(sce))
plot(umap2[ind,],asp=1,col=colorby(sce$condition)[ind], main = 'Condition')
legendby(sce$condition)

# distributions
boxplot(so$AT1_1 ~ sce$clus.epi)



# Proliferation markers
c('Mki67','Top2a')

layout(matrix(1:2,nrow=1))
plot(umap2,asp=1,col=colorby(assay(sce,'binomial_deviance_residuals')['Mki67',], colors = c('grey80','yellow', 'green','blue','darkblue')), main ='Mki67', cex=.5)
plot(umap2,asp=1,col=colorby(assay(sce,'binomial_deviance_residuals')['Top2a',], colors = c('grey80','yellow', 'green','blue','darkblue')), main ='Top2a', cex=.5)
layout(1)


# Cell Cycle? (seems unrelated)
ccgenes <- readRDS('~/Downloads/mouse_cell_cycle_genes_2019.rds')
s.genes <- ccgenes$s.genes
g2m.genes <- ccgenes$m.genes

so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

plot(umap2,asp=1,col=colorby(so$S.Score, colors = c('grey80','yellow', 'green','blue','darkblue')))
plot(umap2,asp=1,col=colorby(so$G2M.Score, colors = c('grey80','yellow', 'green','blue','darkblue')))

# Cell Death? (not seeing much)
apoptosis <- read.delim("~/Downloads/GO_term_apoptotic_process.txt", row.names=NULL)
apop.genes <- apoptosis$MGI.Gene.Marker.ID
apop.epi <- apoptosis[grep('epitheli', apoptosis$Qualifier), ]
apop.epi <- apop.epi[-grep('negative', apop.epi$Qualifier), ]
apop.epi <- apop.epi$MGI.Gene.Marker.ID

so <- AddModuleScore(so, features = list(apop.genes), name = 'apoptosis')
so <- AddModuleScore(so, features = list(apop.epi), name = 'apopEpi')

plot(umap2,asp=1,col=colorby(so$apoptosis1, colors = c('grey80','yellow', 'green','blue','darkblue')))
plot(umap2,asp=1,col=colorby(so$apopEpi1, colors = c('grey80','yellow', 'green','blue','darkblue')))



#####################
# re-add doublets to see if there's a bridge
#####################
# this did not seem to make any difference, just added some more points
require(HDF5Array)
require(SingleCellExperiment)
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
sce <- sce[ ,which(sce$leiden.r1 %in% c(2,14,21,22,32,33,35))]
epi.cells <- colnames(sce)

pca <- readRDS('data/tmp/fastMNNcorrected.rds')
plot(colVars(pca))
dbl <- readRDS('data/tmp/dblclussamp.rds')

# include doublet if any of 5 NNs are in epithelial clus
require(BiocNeighbors)
knn <- findKNN(pca[,1:25], k = 5)
knn <- knn$index[which(dbl$class == 'doublet'),]
keep <- apply(knn,1,function(kn){
  any(rownames(pca)[kn] %in% epi.cells)
})
keep.cells <- rownames(pca)[which(dbl$class == 'doublet')][keep]
cells <- c(epi.cells, keep.cells)

# "re-focus" pca
require(BiocSingular)
pca <- runPCA(pca[cells,], rank=50)
plot(pca$sdev^2)
pca.epi <- pca$x[,1:13]

cd <- colData(sce)
cd$barcode <- rownames(cd)

cd <- data.frame(class = rep('singlet', nrow(pca.epi)))
cd$class[which(cells %in% keep.cells)] <- 'doublet'

require(uwot)
umap2 <- umap(pca.epi, n_components = 2)
umap3 <- umap(pca.epi, n_components = 3)

plot(umap2,asp=1,col=colorby(cd$class), cex=.25)
legendby(cd$class)

plot3d(umap3, aspect = 1, col=colorby(cd$class))

################
################






# RNA velocity for flow between E1/2/3, E7/8
# DE in the Lyz1-high cluster (E4?)
# DE between big AT2 groups (~ RA v. HO)
