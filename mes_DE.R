################
# LOAD / RESET Mesenchymal SCE
################
require(HDF5Array)
require(SingleCellExperiment)
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
anno <- read.csv('data/umap_mes_coordinates.csv', row.names = 1)
anno$mes.clus <- paste0('M',anno$mes.clus)
sce <- sce[, match(rownames(anno), colnames(sce))]
all(colnames(sce) == rownames(anno)) # check
sce$clus.mes <- anno$mes.clus
reducedDim(sce,'umap.mes') <- anno[,c('UMAP1','UMAP2')]
rm(anno)
################
# add logcounts
sce <- logNormCounts(sce)

# get L/R data from CellChat
load('data/CellChatDB.mouse.rda')
db <- CellChatDB.mouse; rm(CellChatDB.mouse)
# db$interaction lists ligand/receptor pairs
# db$geneInfo lists gene Symbol and Synonym (alt. name)

# get list of LIGANDS
ligands <- unique(c(db$interaction$ligand,
                      db$geneInfo$Synonym[db$geneInfo$Symbol %in% db$interaction$ligand],
                      db$geneInfo$Symbol[db$geneInfo$Synonym %in% db$interaction$ligand]))

# Differential Expression
# M1,M3 vs. M2
sub <- sce[, sce$clus.mes %in% c('M1','M2','M3')]

de <- scran::scoreMarkers(sub, groups = sub$clus.mes, block = sub$sample, assay.type = 'logcounts')
de <- de$M2
de$adj.mean.AUC <- de$mean.AUC
de$adj.mean.AUC[de$mean.AUC < .5] <- 1 - de$mean.AUC[de$mean.AUC < .5]
de$ligand <- rownames(de) %in% ligands


# Violin plots
plotExpression(sub, exprs_values = "logcounts",
               features = head(rownames(de[order(de$mean.AUC, decreasing = TRUE), ])), 
               x="condition", colour_by="condition")

plotExpression(sub, exprs_values = "logcounts",
               features = head(rownames(de[order(de$mean.AUC, decreasing = FALSE), ])), 
               x="clus.mes", colour_by="clus.mes")

# ligands only
plotExpression(sub, exprs_values = "logcounts",
               features = head(rownames(de[de$ligand, ][order(de$mean.AUC[de$ligand], decreasing = TRUE), ])), 
               x="clus.mes", colour_by="clus.mes")


# "Volcano" plot
# plot negatives so that directionality is consistent
plot(-de$mean.logFC.cohen, -de$mean.AUC, cex=.5, col=rgb(0,0,0,.1))
# points(-de$mean.logFC.cohen[de$ligand], -de$mean.AUC[de$ligand], cex=.75, col=2)
abline(v=0,lty=2,col='grey'); abline(h=-.5, lty=2,col='grey')
text(-de$mean.logFC.cohen, -de$mean.AUC, labels = rownames(de), cex=.5)






