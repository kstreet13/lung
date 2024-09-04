# create sparse counts matrix from GEO dataset
# counts <- read.csv('~/Downloads/GSE151974_raw_umi_matrix_postfilter.csv.gz')
# rownames(counts) <- counts$X
# counts <- counts[,-1]
# counts <- as.matrix(counts)
# 
# require(Matrix)
# counts <- as(counts,'dgCMatrix')
# saveRDS(counts, file='~/Projects/lung/data/rawcounts.rds')

#counts <- readRDS('~/Projects/lung/data/rawcounts.rds')
#counts <- readRDS('~/OneDrive - University of Southern California/lung_data/rawcounts.rds')
#meta <- read.csv("~/Projects/lung/data/GSE151974_cell_metadata_postfilter.csv.gz")[,-1]
#meta <- read.csv('~/OneDrive - University of Southern California/lung_data/GSE151974_cell_metadata_postfilter.csv.gz')

require(SingleCellExperiment)
#sce <- SingleCellExperiment(assays = list(counts = counts), colData = meta)
#rm(counts, meta)
sce <- readRDS('~/OneDrive - University of Southern California/lung_data/sce.rds')

require(scry)
#sce <- devianceFeatureSelection(sce)

#assay(sce,'logcounts') <- log1p(1000*t(t(assay(sce,'counts')) / colSums(assay(sce,'counts'))))
gene.use <- which(rowData(sce)$binomial_deviance >= sort(rowData(sce)$binomial_deviance, decreasing = TRUE)[2000])
# if deviance doesn't work:
# gene.vars <- rowVars(assay(sce,'logcounts'))
# gene.use <- which(gene.vars >= sort(gene.vars, decreasing = TRUE)[2000])


require(BiocSingular)
#pca <- runPCA(t(assay(sce,'logcounts')[gene.use, ]), rank = 100)
#plot((pca$sdev^2/sum(pca$sdev^2))[1:100])

#reducedDim(sce,'pca') <- pca$x[,1:25]
#rm(pca)

require(uwot)
umap <- umap(reducedDim(sce,'pca'))

plot(umap, asp=1, col = colorby(sce$CellType))
plot(umap, asp=1, col = colorby(assay(sce,'logcounts')[12404, ]))


points(umap[which(sce$CellType=='AT1'), ], col=3)
points(umap[which(sce$CellType=='AT2 1'), ], col=2)
points(umap[which(sce$CellType=='AT2 2'), ], col='firebrick')



subsce <- sce[, which(sce$CellType %in% c('AT1','AT2 1','AT2 2'))]
subsce <- devianceFeatureSelection(subsce)

gene.use <- which(rowData(subsce)$binomial_deviance >= sort(rowData(subsce)$binomial_deviance, decreasing = TRUE)[1000])

subpca <- runPCA(t(assay(subsce,'logcounts')[gene.use, ]), rank = 500)
plot((subpca$sdev^2/sum(subpca$sdev^2))[1:100])

pairs(subpca$x[,1:3], asp=1, col= colorby(subsce$CellType))
pairs(subpca$x[,1:3], asp=1, col= colorby(factor(subsce$Age, levels = c('P3','P7','P14'))))

umap <- umap(subpca$x[,1:20])
plot(umap, asp=1, col = colorby(subsce$CellType))



umap <- umap(pca$x[which(sce$CellType %in% c('AT1','AT2 1','AT2 2','Ciliated','Club')), 1:25])
plot(umap, asp=1, col = colorby(sce$CellType[which(sce$CellType %in% c('AT1','AT2 1','AT2 2','Ciliated','Club'))]))

