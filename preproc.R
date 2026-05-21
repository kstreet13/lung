require(Seurat)
require(SingleCellExperiment)
require(scry)
require(HDF5Array)

s1 <- readRDS('data/Seurat_fromJing/RA_1_filtered_15.rds')
s2 <- readRDS('data/Seurat_fromJing/RA_2_filtered_15.rds')
s3 <- readRDS('data/Seurat_fromJing/RA_3_filtered_15.rds')
s4 <- readRDS('data/Seurat_fromJing/RA_4_filtered_15.rds')
s5 <- readRDS('data/Seurat_fromJing/HO85_1_filtered_15.rds')
s6 <- readRDS('data/Seurat_fromJing/HO85_2_filtered_15.rds')
s7 <- readRDS('data/Seurat_fromJing/HO85_3_filtered_15.rds')
s8 <- readRDS('data/Seurat_fromJing/HO85_4_filtered_15.rds')
# different numbers of genes
ls <- list(s1,s2,s3,s4,s5,s6,s7,s8)
rm(s1,s2,s3,s4,s5,s6,s7,s8)
counts <- lapply(ls, function(so){
  cts <- so@assays$RNA@layers$counts
  rownames(cts) <- rownames(so)
  cn <- gsub('-.*$','',colnames(so))
  cn <- paste0(so$orig.ident[1],'_',cn)
  colnames(cts) <- cn
  return(cts)
})
meta <- lapply(ls, function(so){
  so@meta.data
})
for(i in seq_along(meta)){
  rownames(meta[[i]]) <- colnames(counts[[i]])
}

# get complete list of all genes
genes <- sapply(ls, rownames)
genes <- unique(do.call(c, genes))
rm(ls)

# make all counts matrices the same size
counts <- lapply(counts, function(cts){
  toAdd <- genes[which(! genes %in% rownames(cts))]
  if(length(toAdd) > 0){
    supp <- cts[seq_len(length(toAdd)), ]
    supp[,] <- 0
    rownames(supp) <- toAdd
    cts <- rbind(cts,supp)
  }
  cts <- cts[genes, ]
  return(cts)
})
rm(genes)

# make SCEs
SCEs <- lapply(seq_along(counts), function(i){
  SingleCellExperiment(assay = list(counts = counts[[i]]),
                       colData = meta[[i]])
})
rm(counts, meta)

for(i in seq_along(SCEs)){
  saveHDF5SummarizedExperiment(SCEs[[i]], dir = paste0('data/tmp/sce',i))
}

rm(SCEs)

for(i in 1:8){
  print(i)
  sce <- loadHDF5SummarizedExperiment(paste0('data/tmp/sce',i))
  sce <- nullResiduals(sce, assay="counts", type="deviance")
  saveHDF5SummarizedExperiment(sce, dir = paste0('data/tmp/norm',i))
}

rm(SCEs, sce)


# reset from "norm" objects


sce <- loadHDF5SummarizedExperiment('data/tmp/norm1')
for(i in 2:8){
  sce2 <- loadHDF5SummarizedExperiment(paste0('data/tmp/norm',i))
  sce <- cbind(sce,sce2)
  rm(sce2)
}

# cell-level QC
boxplot(sce$nCount_RNA ~ sce$orig.ident, log = 'y') # HO85_3 looks very different, has some very low-count cells
boxplot(sce$nFeature_RNA ~ sce$orig.ident)

ind <- sample(ncol(sce))
plot(sce$nCount_RNA[ind], sce$nFeature_RNA[ind], cex=.25, col=colorby(sce$orig.ident[ind], alpha=.3))
abline(v = 1600)
abline(h = 800)
# zoom in
plot(sce$nCount_RNA[ind], sce$nFeature_RNA[ind], cex=.25, col=colorby(sce$orig.ident[ind], alpha=.3), xlim = c(0,8000), ylim = c(0,2000))
abline(v = 1600)
abline(h = 800)
# %mito
hist(sce$percent.mt, breaks=200)
abline(v=6)
# define high-quality cells
keep <- which(sce$nCount_RNA > 1600 & sce$nFeature_RNA > 800 & sce$percent.mt < 6)
# percentage of "cells" retained
length(keep) / ncol(sce)
# percentage of counts retained
sum(sce$nCount_RNA[keep]) / sum(sce$nCount_RNA)
# remove low-quality cells
sce <- sce[,keep]

# set up other factors
sce$sample <- sce$orig.ident
sce$condition <- gsub('_.*$', '', sce$orig.ident)

sce <- devianceFeatureSelection(sce, batch = factor(sce$orig.ident), sorted = TRUE)
saveHDF5SummarizedExperiment(sce, dir = 'data/combined8_filt')


# fastMNN - this actually makes a difference
require(batchelor)
sce <- fastMNN(sce, assay.type = 'binomial_deviance_residuals', batch = sce$orig.ident, subset.row = 1:2000) # top 2000 highest deviance genes

saveRDS(reducedDim(sce,'corrected'), file='data/tmp/fastMNNcorrected.rds')

plot(colVars(reducedDim(sce,'corrected')))
require(BiocSingular)
#pca <- runPCA(reducedDim(sce,'corrected'), rank = 50)
#plot(pca$sdev^2)
require(uwot)
umap <- umap(reducedDim(sce,'corrected')[,1:26])

saveRDS(umap, file='data/tmp/fastMNNumap.rds')

# update SCE
rm(sce)
sce <- loadHDF5SummarizedExperiment('data/combined8filt_dimreds/')
reducedDimNames(sce) <- c('OLDpca','OLDumap')
cor <- readRDS('data/tmp/fastMNNcorrected.rds')
reducedDim(sce,'fastMNN') <- cor
rm(cor)
umap <- readRDS('data/tmp/fastMNNumap.rds')
reducedDim(sce,'umap') <- umap
rm(umap)

saveHDF5SummarizedExperiment(sce, dir = 'data/combined8filt_fastMNN')




### OLD STUFF ###



# PCA - probably unnecessary
require(BiocSingular)
pca <- runPCA(t(assay(sce,'binomial_deviance_residuals')[1:2000,]), rank = 50, get.rotation = FALSE)
saveRDS(pca, file = 'data/combined8filt_binomdevresidPCA.rds')

plot(pca$sdev^2)


# UMAP
require(uwot)
pca <- readRDS('data/combined8filt_binomdevresidPCA.rds')
umap <- umap(pca$x[,1:23])

ind <- sample(ncol(sce))
plot(umap[ind,], asp=1, cex=.25, col = colorby(sce$orig.ident[ind]))


# regress out sample?
# l <- lm(pca$x ~ sce$orig.ident)
# l <- l$residuals
# sapply(1:13, function(pc){ cor(l[,pc], pca$x[,pc]) })
# makes no difference


umap <- umap(l[,1:13])

plot(umap[ind,], asp=1, cex=.25, col = colorby(sce$orig.ident[ind]))



reducedDim(sce,'pca') <- pca$x
reducedDim(sce,'umap') <- umap

saveHDF5SummarizedExperiment(sce, dir = 'data/combined8filt_dimreds')

#




