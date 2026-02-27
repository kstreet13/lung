require(SingleCellExperiment)
require(Matrix)

hx <- readRDS('~/OneDrive - University of Southern California/lung_data/rawcounts.rds')
hx <- hx[which(rowSums(hx) > 0), ]
ct <- Seurat::Read10X('~/OneDrive - University of Southern California/lung_data/PN14_Ctrl_filtered_feature_bc_matrix/')
ct <- ct[which(rowSums(ct) > 0), ]

allgenes <- unique(c(rownames(hx),rownames(ct)))
hx.app <- as(matrix(0, nrow = sum(!allgenes %in% rownames(hx)), ncol = ncol(hx)), 'dgCMatrix')
rownames(hx.app) <- allgenes[which(!allgenes %in% rownames(hx))]
hx <- rbind(hx, hx.app)
ct.app <- as(matrix(0, nrow = sum(!allgenes %in% rownames(ct)), ncol = ncol(ct)), 'dgCMatrix')
rownames(ct.app) <- allgenes[which(!allgenes %in% rownames(ct))]
ct <- rbind(ct, ct.app)
rm(hx.app,ct.app)
hx <- hx[allgenes, ]
ct <- ct[allgenes, ]

require(sctransform)
options(future.globals.maxSize = 8000 * 1024^2)
out <- sctransform::vst(hx)

plot(log1p(hx[rownames(out$y),1]), out$y[,1])


require(batchelor)

sce <- batchCorrect(hx, ct, PARAM = FastMnnParam())

reducedDim(sce,'umap') <- uwot::umap(reducedDim(sce,'corrected')[,1:30])

shuf <- sample(ncol(sce))

plot(reducedDim(sce,'umap')[shuf,], asp=1, col=colorby(sce$batch[shuf], alpha=.5), cex=.5)
