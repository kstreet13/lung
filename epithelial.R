# examine epithelial cells, especially AT1/AT2
# 2: AT2
# 14: AT2
# 21: AT1
# 22: AT2
# 32: AT2
# 33: AT1
# 35: Club/Ciliated

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
umaps <- lapply(c(8,14,20,50), function(k){
  umap <- umap(pca$x[,1:k], n_components = 2)
  return(umap)
})
layout(matrix(1:4, 2,2, byrow = TRUE))
for(umap in umaps){
  ind <- sample(nrow(umap))
  plot(umap[ind,],asp=1,col=colorby(as.character(sce$condition[ind])), cex=.25)
  #legendby(as.character(sce$leiden.r1))
}
layout(1)

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

sce$epith.clus <- so$seurat_clusters
#rm(so)

sce$epith.clus <- factor(paste0('E',sce$epith.clus))

plot(umap2,asp=1,col=colorby(as.character(sce$epith.clus)))
legendby(as.character(sce$epith.clus))
labelby(umap2,as.character(sce$epith.clus))


plot3d(umap3, col=colorby(as.character(sce$epith.clus)), aspect = 1)




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
  ClubCil = c('Scgb1a1','Scgb3a2','Foxj1','Dynlrb2')
)
markers$AT1 <- markers$AT1[markers$AT1 %in% rownames(sce)]
markers$AT2 <- markers$AT2[markers$AT2 %in% rownames(sce)]
markers$ClubCil <- markers$ClubCil[markers$ClubCil %in% rownames(sce)]

require(dittoSeq)
dittoDotPlot(sce, assay = 'counts', vars = markers, group.by = 'epith.clus')

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
boxplot(so$AT1_1 ~ sce$epith.clus)



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



# RNA velocity for flow between E1/2/3, E7/8
# re-add doublets to see if there's a bridge

