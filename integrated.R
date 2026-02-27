
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
hx <- SingleCellExperiment(assay = list(counts = hx@assays$RNA@counts[rownames(hx@assays$SCT@scale.data), ], 
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
                    PARAM = FastMnnParam())

reducedDim(sce,'umap') <- uwot::umap(reducedDim(sce,'corrected')[,1:30])

shuf <- sample(ncol(sce))

plot(reducedDim(sce,'umap')[shuf,], asp=1, col=colorby(sce$batch[shuf], alpha=.5), cex=.5)








obj <- merge(hx,ct)
obj <- SCTransform(obj)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
# one-liner to run Integration
obj <- IntegrateLayers(object = obj, method = CCAIntegration,
                       orig.reduction = "pca", new.reduction = 'cca',
                       assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "cca", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "cca_clusters")
obj <- RunUMAP(obj, reduction = "cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("experiment"),
  split.by = 'experiment'
)




obj <- list(hx = hx, ct = ct)

# normalize and identify variable features for each dataset independently
obj <- lapply(X = obj, FUN = function(x) {
  # x <- NormalizeData(x)
  # x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- SCTransform(x, conserve.memory = TRUE)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = obj)

anchors <- FindIntegrationAnchors(object.list = obj, anchor.features = features)

rm(hx,ct,obj)

# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data.
# Note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering
#combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

umap <- combined@reductions$umap@cell.embeddings
lab <- combined@meta.data$orig.ident
lab0 <- lab
lab0[which(lab0 != 'SeuratProject')] <- 'Hyperoxia'
lab0[which(lab0 == 'SeuratProject')] <- 'Control'

shuf <- sample(nrow(umap))

plot(umap[shuf,], asp=1, col = colorby(lab0[shuf], alpha=.5), cex=.5, main='Combined UMAP')
plot.new(); legendby(lab0[shuf], pos='left')

plot(umap, asp=1, col = 'grey90', cex=.5)
points(umap[which(lab0=='Hyperoxia'),], col = brewer.pal(9,'Set1')[2], cex=.5)
points(umap[which(lab0=='Control'),], col = brewer.pal(9,'Set1')[1], cex=.5)

celltype <- combined@meta.data$CellType
plot(umap, asp=1, col = 'grey90', cex=.5)
points(umap, col = c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'))[factor(celltype)], cex=.5)


# inspect area with few Control cells
idx <- which(umap[,1] > -10.5 & umap[,1] < -5.5 &
               umap[,2] > -4 & umap[,2] < 4.5)
plot(umap, asp=1, col = 'grey90', cex=.5, xlim = c(-10.5,-5.5), ylim=c(-4,4.5))
points(umap, col = c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'))[factor(celltype)], cex=.5)
lgnd <- unique(cbind(celltype, c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(12,'Set3'),brewer.pal(9,'Pastel1'))[factor(celltype)]))
lgnd <- lgnd[lgnd[,1] %in% celltype[idx],]
plot.new(); legend('left', col=lgnd[,2], legend=lgnd[,1], pch=16, bty='n')

# grab AT1 and AT2
idx <- which(umap[,1] > -3 & umap[,1] < 7 &
             umap[,2] > 5 & umap[,2] < 11)
idx2 <- which(umap[,1] > -1 &
                umap[,2] >= 11)
ATcells <- c(idx,idx2)
rm(idx,idx2)
ATlab <- celltype[ATcells]
ATlab[is.na(ATlab)] <- 'Control'
ATlab[which(!ATlab %in% c('AT1','AT2 1','AT2 2','Control'))] <- 'Other'

pca <- combined@reductions$pca@cell.embeddings[ATcells,]
pca <- BiocSingular::runPCA(pca, rank=30)


pairs(pca$x[,1:3], asp=1, cex=.5, col = c(3,2,'firebrick',4,'grey')[factor(ATlab)])


subumap <- uwot::umap(pca$x[,1:20])
shuf <- sample(nrow(subumap))
plot(subumap[shuf,], asp=1, cex=.5, col = c(3,2,'firebrick','grey','grey10')[factor(ATlab)][shuf])


