# CellChat
# communication between fibroblasts and epithelial cells
# split by stages of development
# mesenchym ~= myofibroblast
# myofibroblast communication with AT1/AT2
# could compare our P14 to their P14, their P3/P7/P14

#devtools::install_github('jinworks/CellChat')
library(CellChat)

library(SingleCellExperiment)
hx <- readRDS('~/OneDrive - University of Southern California/lung_data/rawcounts.rds')
meta <- read.csv('~/OneDrive - University of Southern California/lung_data/GSE151974_cell_metadata_postfilter.csv.gz')
names(meta)[1] <- 'Cell_ID'
rownames(meta) <- meta$Cell_ID
meta$percent.mt <- NULL # redundant with percent.mito
all(meta$Cell_ID == colnames(hx)) # check
hx <- SingleCellExperiment(assay = list(counts = hx),
                           colData = meta)
rm(meta)

# just one sample
samp <- 'P14_Normoxia'

sce <- hx[, which(hx$Sample == samp)]
assay(sce,'logcounts') <- CellChat::normalizeData(assay(sce,'counts'))
sce$samples <- factor(sce$Sample) # CellChat needs this name, specifically

cch <- createCellChat(object = sce, group.by = "CellType")

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB) # use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# set the used database in the object
cch@DB <- CellChatDB.use

cch <- subsetData(cch) # This step is necessary even if using the whole database (crazy)

cch <- identifyOverExpressedGenes(cch)
cch <- identifyOverExpressedInteractions(cch)


cch <- computeCommunProb(cch, type = "triMean")
cch <- filterCommunication(cch, min.cells = 10)


