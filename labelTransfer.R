# transfer labels from annotated hypoxia data
sce <- readRDS('data/integrated.rds')
hx <- sce[, which(sce$batch == 'hx')]
ct <- sce[, which(sce$batch == 'ct')]

require(BiocNeighbors)
require(SingleCellExperiment)

# find MNNs
mnn <- findMutualNN(reducedDim(ct,'pca'), reducedDim(hx,'pca'), k1 = 30)

# transfer labels
tab <- t(sapply(seq_len(ncol(ct)), function(i1){
  table(hx$CellType[mnn$second[which(mnn$first==i1)]])
}))
ct.celltype <- apply(tab,1,function(x){
  ifelse(max(x)>0, colnames(tab)[which.max(x)], NA)
})
ct.celltype.conf <- apply(tab,1,function(x){
  ifelse(max(x)>0, max(x)/sum(x), 0)
})

res <- data.frame(Cell_ID = ct$Cell_ID,
                  CellType = ct.celltype,
                  CellType_conf = ct.celltype.conf)
res$CellType <- factor(res$CellType, levels = levels(sce$CellType))

x <- merge(colData(sce), res, by = 'Cell_ID',
           all.x = TRUE, sort = FALSE)

