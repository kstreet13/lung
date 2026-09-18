require(HDF5Array)
require(SingleCellExperiment)
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
require(scater)
sce <- logNormCounts(sce)


# get markers (from Jing)
source('celltypemarkers.R')
mdf <- data.frame(gene = unlist(markers), celltype = rep(names(markers), times = lengths(markers)))
rowData(sce)$subMarker <- mdf$celltype[match(rownames(sce), mdf$gene)]
rowData(sce)$Marker <- rowData(sce)$subMarker
rowData(sce)$Marker[rowData(sce)$Marker %in% c('AT1','AT2','Club','Ciliated')] <- 'Epithelial'
rowData(sce)$Marker[rowData(sce)$Marker %in% c('gCAP','aCAP','Artery','Vein','Lymphatic')] <- 'Endothelial'
rowData(sce)$Marker[rowData(sce)$Marker %in% c('Myofibroblast','Pericyte','Smooth muscle','Mesothelial')] <- 'Mesenchymal'

require(dittoSeq)
dittoDotPlot(sce, assay = 'logcounts', vars = markers, group.by = 'leiden.r1')

require(scDotPlot)
scDotPlot(sce, features = unlist(markers), group = 'leiden.r1', featureAnno = 'Marker')

# Annotations of leiden.r1 based on bubble plot
# 1: Immune
# 2: Epithelial - AT2
# 3: Mesenchymal
# 4: Mesenchymal
# 5: Immune
# 6: Endothelial
# 7: Immune
# 8: Immune
# 9: Immune
# 10: Mesenchymal
# 11: Endothelial
# 12: Immune
# 13: Endothelial
# 14: Epithelial - AT2
# 15: Immune
# 16: Mesenchymal
# 17: Endothelial
# 18: Endothelial
# 19: Mesenchymal
# 20: Endothelial
# 21: Epithelial - AT1
# 22: Epithelial - AT2
# 23: Mesenchymal
# 24: Immune
# 25: Immune
# 26: Immune
# 27: Mesenchymal
# 28: Mesenchymal
# 29: ??? Immune?
# 30: Immune
# 31: Mesenchymal
# 32: Epithelial - AT2
# 33: Epithelial - AT1
# 34: Immune
# 35: Epithelial - Club/Ciliated


# UMAP
ind <- sample(ncol(sce))
plot(reducedDim(sce,'umap')[ind,],asp=1,col=colorby(as.character(sce$leiden.r1[ind])), cex=.3)
labelby(reducedDim(sce,'umap'), as.character(sce$leiden.r1))


clus <- sce$leiden.r1

ind <- sample(ncol(sce))
cc <- rep(c(brewer.pal(9,'Set1'), brewer.pal(8,'Set2'), brewer.pal(12,'Set3')), length.out = length(levels(clus)))
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = cc[clus[ind]])

centers <- t(sapply(levels(clus), function(clID){
  colMeans(reducedDim(sce,'umap')[which(clus==clID),])
}))
#legend('right', legend=levels(clus), pch=16, col=cc, bty='n')
points(centers,pch=1,cex=2.5)
points(centers,pch=16,cex=2.5, col=cc)
text(centers, labels = levels(clus), col = 1, font=2)




plot.new()
legend('left',pch=16,col=cc[1:15],legend=levels(sce$leiden.r1)[1:15], bty='n',cex=.5)
legend('center',pch=16,col=cc[16:30],legend=levels(sce$leiden.r1)[16:30], bty='n',cex=.5)
legend('right',pch=16,col=cc[31:35],legend=levels(sce$leiden.r1)[31:35], bty='n',cex=.5)


plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = colorby(sce$condition)[ind])
legendby(sce$condition)

#clus <- 1,14,31
clus <- 31
col <- 'orange'
points(reducedDim(sce,'umap')[which(sce$leiden.r1==clus),], asp=1, cex=.25, col = col, pch=16)

for(clus in c(21,22)){
  col <- 3
  points(reducedDim(sce,'umap')[which(sce$leiden.r1==clus),], asp=1, cex=.25, col = col, pch=16)
}


# by condition (RA/HO)
png("~/Desktop/UMAPbyCondition.png", width = 1000, height = 1000)
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = colorby(as.character(sce$condition[ind])), xlab = 'UMAP-1', ylab = 'UMAP-2')
legendby(as.character(sce$condition), cex=2)
dev.off()

# separated
layout(matrix(1:2,nrow=1))
plot(reducedDim(sce,'umap'), asp=1, cex=.25, col = 'grey90', xlab = 'UMAP-1', ylab = 'UMAP-2', main = 'Room Air')
points(reducedDim(sce,'umap')[which(sce$condition=='RA'),], cex=.25, col = alpha(4,alpha=.3))
plot(reducedDim(sce,'umap'), asp=1, cex=.25, col = 'grey90', xlab = 'UMAP-1', ylab = 'UMAP-2', main = 'Hyperoxia')
points(reducedDim(sce,'umap')[which(sce$condition=='HO85'),], cex=.25, col = alpha(2,alpha=.3))
layout(1)


ind <- sample(ncol(sce))
cc <- rep(c(brewer.pal(9,'Set1'), brewer.pal(8,'Set2'), brewer.pal(12,'Set3')), length.out = length(levels(sce$leiden.r1)))
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = cc[sce$leiden.r1[ind]])


# Marker plot
gene <- 'Mki67'
plot(reducedDim(sce,'umap'), asp=1, cex=.25, col = 'grey80', main=gene)
points(reducedDim(sce,'umap'), asp=1, cex=.25, pch=16,
       col = colorby(assay(sce,'logcounts')[gene,], colors = c('grey80',4,'darkblue')))


# Mesenchymal
plot(reducedDim(sce,'umap'), asp=1, cex=.25, col = 'grey80', main=gene)


# DE to identify ambiguous clusters

require(scran)
de <- findMarkers(sce, groups = sce$leiden.r1, 
                  assay.type = 'binomial_deviance_residuals',
                  test.type = 'wilcox',
                  block = sce$sample)
 


out <- lapply(de, function(x){
  x[1:20,2:4]
})




for(clID in unique(sce$leiden.r1)){
  ind1 <- which(sce$leiden.r1 == clID)
  ind2 <- which(sce$leiden.r1 != clID)
  ct <- assay(sce,'binomial_deviance_residuals')[1,]
  wilcox.test(ct[ind1], ct[ind2])$p.value
  
}



# RA vs HO
cond <- factor(sce$condition, levels = c('RA','HO85'))
tab <- table(cond, sce$leiden.r1)
tab <- t(t(tab) / colSums(tab))
layout(matrix(c(1,2,2),ncol=1))
par(mar=c(.1,4,4,2))
barplot(table(sce$leiden.r1), xlab='')
par(mar=c(5,4,.1,2))
barplot(tab, col = c(4,2))
par(mar=c(5,4,4,2)+.1)
layout(matrix(1))




