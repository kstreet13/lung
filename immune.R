# 1: Immune
# 5: Immune
# 7: Immune
# 8: Immune
# 9: Immune
# 12: Immune
# 15: Immune
# 24: Immune
# 25: Immune
# 26: Immune
# 30: Immune
# 34: Immune

################
# Subset to Immune cells
################
require(HDF5Array)
require(SingleCellExperiment)
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
sce <- sce[, which(sce$leiden.r1 %in% c(1,5,7,8,9,12,15,24,25,26,30,34))]
################


require(BiocSingular)
pca.imm <- runPCA(reducedDim(sce,'fastMNN'), rank = 50)
plot(pca.imm$sdev^2) # non-linear: ~16

pairs(pca.imm$x[,1:4], col=colorby(as.character(sce$leiden.r1)), asp=1, cex=.25)

options(rgl.useNULL = TRUE)
library(rgl)
options(rgl.printRglwidget = TRUE)

plot3d(pca.imm$x[,1:3], col=colorby(as.character(sce$leiden.r1)), aspect = 1)


require(uwot)
umap2 <- umap(pca.imm$x[,1:16], n_components = 2)
umap3 <- umap(pca.imm$x[,1:16], n_components = 3)

plot(umap2,asp=1,col=colorby(as.character(sce$leiden.r1)))
legendby(as.character(sce$leiden.r1))

plot3d(umap3, col=colorby(as.character(sce$leiden.r1)), aspect = 1)


##############
# RE-CLUSTER #
##############
# perform Leiden clustering in Seurat
# set up Seurat object
require(Seurat)
require(Matrix)
so <- Matrix(0, nrow = nrow(sce), ncol = ncol(sce), sparse = TRUE)
so <- CreateSeuratObject(so)
pca <- CreateDimReducObject(embeddings = pca.imm$x[,1:16], key = "PC_")
colnames(so) <- rownames(pca)
so@reductions[['pca']] <- pca
rm(pca)

# clustering
so <- FindNeighbors(so, reduction = 'pca')
so <- FindClusters(so, algorithm = 4, resolution = .3)

sce$clus.imm <- so$seurat_clusters
#rm(so)
sce$clus.imm <- factor(paste0('imm',sce$clus.imm))
levels(sce$clus.imm) <- paste0('imm',1:length(unique(sce$clus.imm)))

# plot
plot(umap2,asp=1,col=colorby(as.character(sce$clus.imm)))
#legendby(as.character(sce$clus.imm))
labelby(umap2,as.character(sce$clus.imm))

################
# SAVE IMMUNE STUFF (DR/CLUS)
################
saveRDS(list(clus.imm = sce$clus.imm,
             pca.imm = pca.imm$x[,1:16],
             umap2.imm = umap2,
             umap3.imm = umap3),
        file = 'data/immuneANNO.rds')
###

################
# LOAD / RESET Immune SCE
################
require(HDF5Array)
require(SingleCellExperiment)
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
anno <- readRDS('data/immuneANNO.rds')
sce <- sce[, match(rownames(anno$pca.imm), colnames(sce))]
sce$clus.imm <- anno$clus.imm
reducedDim(sce,'pca.imm') <- anno$pca.imm
reducedDim(sce,'umap.imm') <- anno$umap2.imm
reducedDim(sce,'umap3.imm') <- anno$umap3.imm
rm(anno)
################
################

# UMAP plots
ind <- sample(ncol(sce))
plot(reducedDim(sce,'umap.imm')[ind,],asp=1,col=colorby(as.character(sce$clus.imm[ind])), cex=.5)
labelby(reducedDim(sce,'umap.imm'),as.character(sce$clus.imm))

plot(reducedDim(sce,'umap.imm')[ind,],asp=1,col=colorby(sce$condition)[ind], cex=.5)
legendby(sce$condition)

layout(matrix(1:2,nrow=1))
ind <- which(sce$condition=='RA')
plot(reducedDim(sce,'umap.imm')[ind,],asp=1,col=colorby(sce$condition)[ind], cex=.25)
ind <- which(sce$condition=='HO85')
plot(reducedDim(sce,'umap.imm')[ind,],asp=1,col=colorby(sce$condition)[ind], cex=.25)


# cluster x condition
barplot(table(sce$condition,sce$clus.imm), col=c(2,4), las=2, beside=TRUE)
barplot(table(sce$condition,sce$clus.imm), col=c(2,4), las=2, beside=FALSE)

layout(1)

# marker plots
gene <- 'Mki67'
plot(reducedDim(sce,'umap.imm')[ind,], asp=1,
     col=colorby(assay(sce,'binomial_deviance_residuals')[gene,ind],
                 colors = c('grey90','grey90','lightgreen','green','darkgreen','blue','darkblue')),
     cex=.5, main=gene)


cv <- sapply(sortu(sce$clus.imm), function(cl){ colorby(sce$clus.imm)[which.max(sce$clus.imm==cl)] })

boxplot(assay(sce,'binomial_deviance_residuals')[gene,] ~ sce$clus.imm)

source('~/Projects/OLD/thingsandstuff/violinplot.R')
violinplot(by(assay(sce,'binomial_deviance_residuals')[gene,], sce$clus.imm, c), las=2, col=cv)


############################
# CELL TYPE IDENTIFICATION #
############################

# Marker genes
# based on: https://www.cellsignal.com/pathways/immune-cell-markers-mouse
{
  # Marker names were converted from proteins to mouse gene symbols:
  #   CD45 -> Ptprc
  #   Cd11b -> Itgam
  #   F4/80 -> Adgre1
  #   Ly-6G -> Ly6g
  #   CD317 -> Bst2
  #   NKp46 -> Ncr1
  #   NK1.1 -> Klrb1c
  #   Cd138 -> Sdc1
  #   Bcma -> Tnfrsf17
  #   Cd11c -> Itgax
  # - The guide was designed for IHC, so these markers should be treated as
  #   candidate scRNA-seq markers rather than definitive annotation rules.
  # - Negative markers can be unreliable in scRNA-seq because of dropout.
  # - Klrb1c/NK1.1 expression is mouse-strain dependent.
  marker_genes <- data.frame(
    cell_class = c(
      "Leukocyte",
      
      rep("Myeloid", 35),
      
      rep("Lymphoid", 49)
    ),
    
    cell_type = c(
      "Leukocyte",
      
      "Myeloid cell",
      "Neutrophil", "Neutrophil",
      "Eosinophil", "Eosinophil", "Eosinophil",
      "Basophil", "Basophil",
      "Mast cell", "Mast cell", "Mast cell",
      "Monocytic MDSC", "Monocytic MDSC", "Monocytic MDSC",
      "Polymorphonuclear MDSC", "Polymorphonuclear MDSC",
      "Polymorphonuclear MDSC", "Polymorphonuclear MDSC",
      "Plasmacytoid dendritic cell", "Plasmacytoid dendritic cell",
      "Activated pDC",
      "Conventional dendritic cell", "Conventional dendritic cell",
      "Activated cDC", "MregDC",
      "Langerhans cell",
      "cDC1", "cDC1",
      "Resident cDC1",
      "Migratory cDC1",
      "cDC2", "cDC2",
      "Macrophage",
      "M1-like macrophage", "M1-like macrophage",
      
      "M1-like macrophage",
      "M2-like macrophage", "M2-like macrophage", "M2-like macrophage",
      "Monocyte", "Monocyte",
      "Alpha-beta T cell",
      "Cytotoxic T cell",
      "Helper T cell",
      "Th1", "Th1",
      "Th2", "Th2",
      "Th9", "Th9",
      "Tfh", "Tfh", "Tfh",
      "Th17", "Th17",
      "Th22", "Th22",
      "Regulatory T cell", "Regulatory T cell",
      "Gamma-delta T cell", "Gamma-delta T cell",
      "Gamma-delta T cell", "Gamma-delta T cell",
      "NKT cell", "NKT cell",
      "Type I NKT cell", "Type I NKT cell",
      "NK cell", "NK cell", "NK cell",
      "NK cell", "NK cell",
      "Activated NK cell",
      "Cytotoxic NK cell", "Cytotoxic NK cell",
      "B cell",
      "Naive B cell", "Naive B cell",
      "Switched-memory B cell", "Switched-memory B cell",
      "Unswitched-memory B cell", "Unswitched-memory B cell",
      "Plasma cell", "Plasma cell"
    ),
    
    gene = c(
      "Ptprc",
      
      "Itgam",
      "Adgre1", "Ly6g",
      "Ccr3", "Adgre1", "Siglecf",
      "Fcer1a", "Kit",
      "Fcer1a", "Kit", "Fcer2a",
      "Ly6c2", "Ly6g", "Arg1",
      "Ly6c2", "Ly6g", "Arg1", "Adgre1",
      "Siglech", "Bst2",
      "Cd83",
      "Itgax", "H2-Ab1",
      "Cd83", "Lamp3",
      "Cd207",
      "Xcr1", "Clec9a",
      "Cd8a",
      "Itgae",
      "Itgam", "Sirpa",
      "Adgre1",
      "Cd86", "Cd80",
      
      "Nos2",
      "Cd163", "Mrc1", "Arg1",
      "Cd14", "Adgre1",
      "Cd3e",
      "Cd8a",
      "Cd4",
      "Tbx21", "Ifng",
      "Gata3", "Il4",
      "Spi1", "Il9",
      "Bcl6", "Cxcr5", "Il21",
      "Rorc", "Il17a",
      "Ahr", "Il22",
      "Foxp3", "Il2ra",
      "Cd4", "Cd8a", "Cd3e", "Trdc",
      "Cd3e", "Klrb1c",
      "Trav11", "Traj18",
      "Klrb1c", "Ncr1", "Klrk1",
      "Cd3e", "Klrb1c",
      "Cd69",
      "Gzmb", "Prf1",
      "Cd19",
      "Ighd", "Cd27",
      "Ighd", "Cd27",
      "Ighd", "Cd27",
      "Tnfrsf17", "Sdc1"
    ),
    
    expression = c(
      "positive",
      
      "positive",
      "medium", "positive",
      "positive", "medium", "positive",
      "positive", "negative",
      "positive", "positive", "positive",
      "positive", "negative", "positive",
      "low", "positive", "positive", "not specified",
      "positive", "positive",
      "positive",
      "positive", "positive",
      "positive", "positive",
      "positive",
      "positive", "positive",
      "positive",
      "positive",
      "positive", "positive",
      "high",
      "positive", "positive",
      
      "positive",
      "positive", "positive", "positive",
      "positive", "negative/low",
      "positive",
      "positive",
      "positive",
      "positive", "positive",
      "positive", "positive",
      "positive", "positive",
      "positive", "positive", "positive",
      "positive", "positive",
      "positive", "positive",
      "positive", "positive",
      "negative", "negative", "positive", "positive",
      "positive", "positive",
      "positive", "positive",
      "positive", "positive", "positive",
      "negative", "positive",
      "positive",
      "positive", "positive",
      "positive",
      "positive", "negative",
      "negative", "positive",
      "positive", "positive",
      "positive", "positive"
    ),
    stringsAsFactors = FALSE
  )
  marker_genes <- marker_genes[order(marker_genes$cell_class, marker_genes$cell_type, marker_genes$gene), ]
  row.names(marker_genes) <- NULL
}
head(marker_genes)
marker_genes <- marker_genes[marker_genes$gene %in% rownames(sce), ]
marker_genes$expressed <- sapply(marker_genes$gene, function(g){
  sum(assay(sce,'counts')[g,] > 0) >= 20
})


sce <- logNormCounts(sce)


rowData(sce)$cellclassMarker <- marker_genes$cell_class[match(rownames(sce), marker_genes$gene)]
rowData(sce)$celltypeMarker <- marker_genes$cell_type[match(rownames(sce), marker_genes$gene)]
rowData(sce)$directionMarker <- marker_genes$expression[match(rownames(sce), marker_genes$gene)]

toplevel <- c('Ptprc','Itgam','Adgre1','Fcer1a','Arg1','Bst2','Siglech','Itgax','H2-Ab1','Cd207','Cd14','Cd3e','Cd4','Cd8a','Trdc','Klrb1c','Klrk1','Ncr1','Cd19','Sdc1','Tnfrsf17','Cd3d','Nkg7','Cd79a','Cd68','Gzmb','Cd3g','Cd79b')


require(scDotPlot)
png("~/Desktop/dots.png", width = 900, height = 2000, res=200)
scDotPlot(sce, features = unique(marker_genes$gene[marker_genes$expressed]), group = 'clus.imm')
dev.off()

require(scDotPlot)
png("~/Desktop/dots.png", width = 900, height = 1500, res=200)
scDotPlot(sce, features = toplevel, group = 'clus.imm')
dev.off()






