require(HDF5Array)
require(SingleCellExperiment)

sce <- loadHDF5SummarizedExperiment('data/combined8filt_leiden/')

# perform Leiden clustering in Seurat
# set up Seurat object
require(Seurat)
require(Matrix)
so <- Matrix(0, nrow = nrow(sce), ncol = ncol(sce), sparse = TRUE)
so <- CreateSeuratObject(so)
pca <- CreateDimReducObject(embeddings = reducedDim(sce,'fastMNN'), key = "PC_")
colnames(so) <- rownames(pca)
so@reductions[['pca']] <- pca
rm(pca)

# clustering
so <- FindNeighbors(so, reduction = 'pca')
so <- FindClusters(so, algorithm = 4, resolution = 1)

sce$leiden.r1 <- so$seurat_clusters
#rm(so)

ind <- sample(ncol(sce))
cc <- rep(c(brewer.pal(9,'Set1'), brewer.pal(8,'Set2'), brewer.pal(12,'Set3')), length.out = length(levels(sce$leiden.r1)))
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = cc[sce$leiden.r1[ind]])


# get markers (from Jing)
source('celltypemarkers.R')

# try to find appropriate ordering for clusters
sce$clus <- sce$leiden.r1
means <- sapply(levels(sce$clus), function(clID){
  colMeans(reducedDim(sce,'fastMNN')[which(sce$clus==clID), ])
})
ord <- hclust(dist(t(means)))$order

levels(sce$clus) <- ord
stopifnot(all(!is.na(sce$clus)))

require(dittoSeq)
dittoDotPlot(sce, assay = 'counts', vars = markers, group.by = 'clus')


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

