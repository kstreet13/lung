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
plot(umap, asp=1, col = c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'))[factor(sce$CellType)])
plot(umap, asp=1, col = colorby(assay(sce,'logcounts')[12404, ])) #sftpc

plot(umap, asp=1, col = colorby(assay(sce,'logcounts')['Nkx2-1', ]), main='Nkx2.1') #epithelial
plot(umap, asp=1, col = colorby(assay(sce,'logcounts')['Epcam', ]), main='Epcam') #epithelial

points(umap[which(sce$CellType=='AT1'), ], col=3)
points(umap[which(sce$CellType=='AT2 1'), ], col=2)
points(umap[which(sce$CellType=='AT2 2'), ], col='firebrick')


# just AT1 and AT2 cells
#########################
subsce <- sce[, which(sce$CellType %in% c('AT1','AT2 1','AT2 2'))]

subsce <- devianceFeatureSelection(subsce)
gene.use <- which(rowData(subsce)$binomial_deviance >= sort(rowData(subsce)$binomial_deviance, decreasing = TRUE)[1000])

subpca <- runPCA(t(assay(subsce,'logcounts')[gene.use, ]), rank = 50)
plot((subpca$sdev^2/sum(subpca$sdev^2))[1:50])

pairs(subpca$x[,1:3], asp=1, col= colorby(subsce$CellType))
pairs(subpca$x[,1:3], asp=1, col= colorby(factor(subsce$Age, levels = c('P3','P7','P14'))))

subumap <- umap(subpca$x[,1:20])
plot(subumap, asp=1, col = colorby(subsce$CellType))

plot(subumap, asp=1, col = 'grey80')
points(subumap[which(assay(subsce,'counts')['Mki67',] > 0), ], col = 2)
points(subumap[which(assay(subsce,'counts')['Top2a',] > 0), ], col = 2)



#umap <- umap(pca$x[which(sce$CellType %in% c('AT1','AT2 1','AT2 2','Ciliated','Club')), 1:25])
#plot(umap, asp=1, col = colorby(sce$CellType[which(sce$CellType %in% c('AT1','AT2 1','AT2 2','Ciliated','Club'))]))





rd <- subpca$x[subsce$CellType %in% c('AT2 1','AT2 2'), 1:3]
cl <- subsce$CellType[subsce$CellType %in% c('AT2 1','AT2 2')]

require(slingshot)
pst <- slingshot(rd)

plot(subpca$x[,1:2], asp=1, col= colorby(subsce$CellType))
points(subpca$x[subsce$CellType %in% c('AT2 1','AT2 2'),1:2], col = colorby(slingPseudotime(pst)[,1]))

plot(subumap, asp=1, col = colorby(subsce$CellType))
points(subumap[subsce$CellType %in% c('AT2 1','AT2 2'), ], col = colorby(slingPseudotime(pst)[,1]))



expr <- assay(subsce,'logcounts')[,subsce$CellType %in% c('AT2 1','AT2 2')]
expr <- expr[rowVars(expr) > 0, ]

pstCors <- apply(expr,1,function(g){
  cor(g, slingPseudotime(pst)[,1], method='spearman')
})


which.max(abs(pstCors))

plot(slingPseudotime(pst)[,1], assay(subsce,'logcounts')['Rpl23',subsce$CellType %in% c('AT2 1','AT2 2')])

plot(subumap, asp=1, col = colorby(assay(subsce,'logcounts')['Rpl23',]))


pstCors[order(abs(pstCors), decreasing = TRUE)]



