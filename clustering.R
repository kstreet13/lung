require(HDF5Array)
require(SingleCellExperiment)

sce <- loadHDF5SummarizedExperiment('data/combined8filt_dbl/')
sce <- sce[ ,which(sce$dbl.clus.samp$class == 'singlet')]
sce$leiden.r1.predbl <- sce$leiden.r1

# re-run PCA and UMAP after removing doublets
sce <- loadHDF5SummarizedExperiment('data/combined8filt/')
reducedDim(sce,'OLDpca') <- NULL
reducedDim(sce,'OLDumap') <- NULL
reducedDim(sce,'umap_dbl') <- reducedDim(sce,'umap')
reducedDim(sce,'umap') <- NULL

require(BiocSingular)
pca <- runPCA(reducedDim(sce,'fastMNN'), rank = 50)
#plot(pca$sdev^2)
reducedDim(sce,'pca') <- pca$x
rm(pca)
require(uwot)
reducedDim(sce,'umap') <- umap(reducedDim(sce,'pca')[,1:24])

# perform Leiden clustering in Seurat (somehow, clusters are identical whether you use 24 or 50 PCs)
# set up Seurat object
require(Seurat)
require(Matrix)
so <- Matrix(0, nrow = nrow(sce), ncol = ncol(sce), sparse = TRUE)
so <- CreateSeuratObject(so)
pca <- CreateDimReducObject(embeddings = reducedDim(sce,'pca'), key = "PC_")
colnames(so) <- rownames(pca)
so@reductions[['pca']] <- pca
rm(pca)

# clustering
so <- FindNeighbors(so, reduction = 'pca')
so <- FindClusters(so, algorithm = 4, resolution = 1)

sce$leiden.r1.oldpca <- sce$leiden.r1
sce$leiden.r1 <- so$seurat_clusters
#rm(so)

ind <- sample(ncol(sce))
cc <- rep(c(brewer.pal(9,'Set1'), brewer.pal(8,'Set2'), brewer.pal(12,'Set3')), length.out = length(levels(sce$leiden.r1)))
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = cc[sce$leiden.r1[ind]])

# save
saveHDF5SummarizedExperiment(sce, dir='data/combined8reclus')


# get markers (from Jing)
source('celltypemarkers.R')

sce <- loadHDF5SummarizedExperiment('data/combined8reclus')

require(dittoSeq)
dittoDotPlot(sce, assay = 'counts', vars = markers, group.by = 'leiden.r1')


# doublets?
boxplot(sce$dbl.noclus.samp$score ~ sce$clus)
mosaicplot(table(sce$clus, sce$dbl.noclus.samp$class), col=2:1)


ind <- sample(ncol(sce))

# UMAP plots grouped by sample
png("~/Desktop/UMAPbySample.png", width = 1000, height = 1000)
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = colorby(as.character(sce$sample[ind])), xlab = 'UMAP-1', ylab = 'UMAP-2')
legendby(as.character(sce$sample))
dev.off()

# condition
png("~/Desktop/UMAPbyCondition.png", width = 1000, height = 1000)
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = colorby(as.character(sce$condition[ind])), xlab = 'UMAP-1', ylab = 'UMAP-2')
legendby(as.character(sce$condition))
dev.off()

# clustering
png("~/Desktop/UMAPbyCluster.png", width = 1000, height = 1000)
cc <- rep(c(brewer.pal(9,'Set1'), brewer.pal(8,'Set2'), brewer.pal(12,'Set3')), length.out = length(levels(sce$clus)))
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = cc[sce$clus[ind]], xlab = 'UMAP-1', ylab = 'UMAP-2')
dev.off()

# UMAP grouped by sample that shows only the RA samples
png("~/Desktop/UMAPbySampleRAonly.png", width = 1000, height = 1000)
plot(reducedDim(sce,'umap'), asp=1, cex=.25, col = 'grey80', xlab = 'UMAP-1', ylab = 'UMAP-2')
raind <- sample(which(sce$condition == 'RA'))
points(reducedDim(sce,'umap')[raind,], cex=.25, col = colorby(as.character(sce$sample))[raind])
legendby(as.character(sce$sample))
dev.off()

