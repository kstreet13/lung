# transfer labels from annotated hypoxia data
require(SingleCellExperiment)

# preprocessing "their" dataset (from Nature paper?)
theirs <- readRDS('~/OneDrive - University of Southern California/lung_data/rawcounts.rds')
theirs <- SingleCellExperiment(assays = list(counts = theirs))
int <- readRDS('data/integrated.rds')
stopifnot(all(colnames(theirs) %in% int$Cell_ID))
colData(theirs) <- colData(int)[match(colnames(theirs), int$Cell_ID), ]
rm(int)
theirs$Sample <- factor(as.character(theirs$Sample))
require(scry)
theirs <- devianceFeatureSelection(theirs, batch = theirs$Sample, sorted = TRUE)
#theirs <- nullResiduals(theirs, batch = theirs$Sample)
l <- lapply(levels(theirs$Sample), function(samp){
  sub <- theirs[, which(theirs$Sample == samp)]
  sub <- nullResiduals(sub)
  return(sub)
})
rm(theirs)
#x <- do.call(cbind, l) # crashes
theirs <- l[[1]]
l <- l[-1]
ll <- length(l)
for(ii in 1:ll){
  theirs <- cbind(theirs, l[[1]])
  l <- l[-1]
  print(length(l))
}
rm(l,ll,ii)
top2Kish <- rownames(theirs)[1:2000]
top2Ktheirs <- top2Kish

# load our dataset
require(HDF5Array)
sce <- loadHDF5SummarizedExperiment('data/combined8filt_fastMNN/')

# get top 2K genes from both, use union
top2Kish <- unique(c(top2Kish, rownames(sce)[1:2000]))
top2Kish <- top2Kish[top2Kish %in% rownames(sce) &
                       top2Kish %in% rownames(theirs)]

theirs <- theirs[top2Ktheirs, ]

require(SingleR)
pred.hesc <- SingleR(test = sce, ref = theirs, labels = theirs$CellType,
                     assay.type.test=2,
                     assay.type.ref=2)








# x = SCE built with just binom_dev_resids from 'sce'
# y = SCE built with just binom_dev_resids from 'theirs'
# merge those (because it seems to work with SCE objects, but not matrices)

fullmerge <- cbind(assay(sce,'binomial_deviance_residuals')[top2Kish,],
                   assay(theirs,'binomial_deviance_residuals')[top2Kish,])




require(BiocNeighbors)

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
                  infCellType = ct.celltype,
                  infCellType_conf = ct.celltype.conf)
res$infCellType <- factor(res$infCellType, levels = levels(sce$CellType))

x <- merge(colData(sce), res, by = 'Cell_ID',
           all.x = TRUE, sort = FALSE)
x$combinedCellType <- ifelse(is.na(x$CellType), x$infCellType, x$CellType)
x <- x[match(colnames(sce), x$Cell_ID), ]

# check
all(colnames(sce) == x$Cell_ID)

colData(sce) <- x

saveRDS(sce, file = 'data/integrated.rds')








# Evaluation

# Epithelial cells: Nkx2-1
# Endothelial cells: Cdh5
# Mesenchymal cells: Col1a1
# Immune cells: Ptprc
# 
# Epithelial subtypes
# AT1: Hopx
# AT2: Sftpc
# Club: Scgb1a1
# Airway: Foxj1
# 
# Mesenchyme subtypes
# Myofibroblast: Pdgfra Inmt Fn1
# Pericyte: Pdgfrb Ebf1 Postn
# Smooth muscle: Tagln Acta2

markers <- c('Nkx2-1','Cdh5','Col1a1','Ptprc','Hopx','Sftpc','Scgb1a1','Foxj1',
             'Pdgfra','Inmt','Fn1','Pdgfrb','Ebf1','Postn','Tagln','Acta2')

require(scuttle)
require(dittoSeq)

hx <- readRDS('~/OneDrive - University of Southern California/lung_data/rawcounts.rds')
meta <- read.csv('~/OneDrive - University of Southern California/lung_data/GSE151974_cell_metadata_postfilter.csv.gz')
names(meta)[1] <- 'Cell_ID'
rownames(meta) <- meta$Cell_ID
all(meta$Cell_ID == colnames(hx)) # check
hx <- SingleCellExperiment(assay = list(counts = hx),
                           colData = meta)
rm(meta)
hx <- logNormCounts(hx)

ct <- Seurat::Read10X('~/OneDrive - University of Southern California/lung_data/PN14_Ctrl_filtered_feature_bc_matrix/')
ct <- SingleCellExperiment(assays = list(counts = ct))
ct <- logNormCounts(ct)
all(colnames(ct) == sce$Cell_ID[which(sce$Cell_ID %in% colnames(ct))]) # check
ct$infCellType <- sce$infCellType[which(sce$Cell_ID %in% colnames(ct))]
ct$infCellType <- as.character(ct$infCellType)
ct$infCellType[is.na(ct$infCellType)] <- 'NA'
ct$infCellType <- factor(ct$infCellType)

dittoDotPlot(hx, vars = markers[markers%in%rownames(hx)], group.by = 'CellType', vars.dir = 'x')

dittoDotPlot(ct, vars = markers[markers%in%rownames(ct)], group.by = 'infCellType', vars.dir = 'x')
