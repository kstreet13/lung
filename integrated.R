
# compare their data with "our" data
# PN14 Ctrl (ours) vs. hyperoxia (theirs, sce from GEO)

require(Seurat)
options(future.globals.maxSize = 16000 * 1024^2)
require(SingleCellExperiment)

# get appropriate set of HVGs
{
  # hx <- readRDS('~/OneDrive - University of Southern California/lung_data/rawcounts.rds')
  # hx <- CreateSeuratObject(counts = hx)
  # hx <- UpdateSeuratObject(hx)
  # hx <- SCTransform(hx)
  # # extract genes
  # hvg.hx <- rownames(hx@assays$SCT@scale.data)
  # rm(hx)
  # 
  # ct <- Read10X('~/OneDrive - University of Southern California/lung_data/PN14_Ctrl_filtered_feature_bc_matrix/')
  # ct <- CreateSeuratObject(counts = ct)
  # ct <- UpdateSeuratObject(ct)
  # ct <- SCTransform(ct)
  # # extract genes
  # hvg.ct <- rownames(ct@assays$SCT@scale.data)
  # rm(ct)
  # 
  # hvg <- unique(c(hvg.hx, hvg.ct))
  # saveRDS(hvg, file='data/hvg_for_integration.rds')
  # rm(hvg.hx,hvg.ct)
}

hvgs <- readRDS('data/hvg_for_integration.rds')

hx <- readRDS('~/OneDrive - University of Southern California/lung_data/rawcounts.rds')
hx <- CreateSeuratObject(counts = hx)
hx <- UpdateSeuratObject(hx)
hx <- SCTransform(hx, residual.features = hvgs, variable.features.n = NULL)
hx <- SingleCellExperiment(assay = list(counts = hx@assays$SCT@counts[rownames(hx@assays$SCT@scale.data), ], 
                                        resids = hx@assays$SCT@scale.data))

ct <- Read10X('~/OneDrive - University of Southern California/lung_data/PN14_Ctrl_filtered_feature_bc_matrix/')
ct <- CreateSeuratObject(counts = ct)
ct <- UpdateSeuratObject(ct)
ct <- SCTransform(ct, residual.features = hvgs, variable.features.n = NULL)
ct <- SingleCellExperiment(assay = list(counts = ct@assays$SCT@counts[rownames(ct@assays$SCT@scale.data), ], 
                                        resids = ct@assays$SCT@scale.data))

require(batchelor)
common.genes <- rownames(hx)[rownames(hx) %in% rownames(ct)]

sce <- batchCorrect(hx[common.genes, ], ct[common.genes],
                    assay.type = 'resids',
                    PARAM = ClassicMnnParam())
sce$batch <- factor(sce$batch)
levels(sce$batch) <- c('hx','ct')

reducedDim(sce,'pca') <- BiocSingular::runPCA(t(assay(sce,'corrected')), rank = 30)$x
reducedDim(sce,'umap') <- uwot::umap(reducedDim(sce,'pca')[,1:30])

saveRDS(sce, file = 'data/integrated.rds')


shuf <- sample(ncol(sce))

plot(reducedDim(sce,'umap')[shuf,], asp=1, col=colorby(sce$batch[shuf], alpha=.5), cex=.5)

layout(matrix(1:2, nrow=1))
ind <- which(sce$batch=='hx')
plot(reducedDim(sce,'umap'), asp=1, col='grey90', cex=.5, main='Hypoxia')
points(reducedDim(sce,'umap')[ind,], col=colorby(sce$batch, alpha=.5)[ind], cex=.5)
ind <- which(sce$batch=='ct')
plot(reducedDim(sce,'umap'), asp=1, col='grey90', cex=.5, main='Control (Ours)')
points(reducedDim(sce,'umap')[ind,], col=colorby(sce$batch, alpha=.5)[ind], cex=.5)


# add metadata
sce <- readRDS('data/integrated.rds')
colData(sce)$Cell_ID <- colnames(sce)
# load existing metadata
meta <- read.csv('~/OneDrive - University of Southern California/lung_data/GSE151974_cell_metadata_postfilter.csv.gz')
names(meta)[1] <- 'Cell_ID'
meta$percent.mt <- NULL # redundant with percent.mito
# add easy stuff from our data
ct <- Read10X('~/OneDrive - University of Southern California/lung_data/PN14_Ctrl_filtered_feature_bc_matrix/')
ct <- CreateSeuratObject(counts = ct)
ct$Cell_ID <- colnames(ct)
ct$Age <- 'P14'
ct$Oxygen <- 'Normoxia'
ct$Sample <- 'P14_Normoxia_Ours'
meta <- merge(meta, ct@meta.data, all.x = TRUE, all.y = TRUE, sort = FALSE)

# check before merging
all(meta$Cell_ID %in% colnames(sce))
all(colnames(sce) %in% meta$Cell_ID)
# merge (preserving colnames(sce) and rownames(colData(sce)))
cd <- merge(colData(sce), meta, by = 'Cell_ID', all.x = TRUE, sort = FALSE)
rownames(cd) <- cd$Cell_ID
colData(sce) <- cd
rm(ct,cd,meta)
# convert to factors
sce$orig.ident <- factor(sce$orig.ident)
sce$Age <- factor(sce$Age)
sce$Oxygen <- factor(sce$Oxygen)
sce$Phase <- factor(sce$Phase)
sce$Sample <- factor(sce$Sample)
sce$cluster_high_res <- factor(sce$cluster_high_res)
sce$CellType <- factor(sce$CellType)

saveRDS(sce, file = 'data/integrated.rds')



### OLD STUFF ###

# inspect area with few Control cells
# idx <- which(umap[,1] > -10.5 & umap[,1] < -5.5 &
#                umap[,2] > -4 & umap[,2] < 4.5)
# plot(umap, asp=1, col = 'grey90', cex=.5, xlim = c(-10.5,-5.5), ylim=c(-4,4.5))
# points(umap, col = c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'))[factor(celltype)], cex=.5)
# lgnd <- unique(cbind(celltype, c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'))[factor(celltype)]))
# lgnd <- lgnd[lgnd[,1] %in% celltype[idx],]
# plot.new(); legend('left', col=lgnd[,2], legend=lgnd[,1], pch=16, bty='n')

# grab AT1 and AT2
# idx <- which(umap[,1] > -3 & umap[,1] < 7 &
#              umap[,2] > 5 & umap[,2] < 11)
# idx2 <- which(umap[,1] > -1 &
#                 umap[,2] >= 11)
# ATcells <- c(idx,idx2)
# rm(idx,idx2)
# ATlab <- celltype[ATcells]
# ATlab[is.na(ATlab)] <- 'Control'
# ATlab[which(!ATlab %in% c('AT1','AT2 1','AT2 2','Control'))] <- 'Other'

# pca <- combined@reductions$pca@cell.embeddings[ATcells,]
# pca <- BiocSingular::runPCA(pca, rank=30)


# pairs(pca$x[,1:3], asp=1, cex=.5, col = c(3,2,'firebrick',4,'grey')[factor(ATlab)])


# subumap <- uwot::umap(pca$x[,1:20])
# shuf <- sample(nrow(subumap))
# plot(subumap[shuf,], asp=1, cex=.5, col = c(3,2,'firebrick','grey','grey10')[factor(ATlab)][shuf])


