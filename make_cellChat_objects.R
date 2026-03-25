# make CellChat objects that will be useful

# 1. Our data (P14)
# 2. Their data - P14 Normoxia
# 3. Their data - P7 Normoxia
# 4. Their data - P3 Normoxia
library(CellChat)
library(SingleCellExperiment)

### OUR DATA
ct <- Seurat::Read10X('~/OneDrive - University of Southern California/lung_data/PN14_Ctrl_filtered_feature_bc_matrix/')
int <- readRDS('data/integrated.rds')
int <- int[, which(int$Cell_ID %in% colnames(ct))]
all(int$Cell_ID == colnames(ct)) # check
sce <- SingleCellExperiment(assays = list(counts = ct),
                           colData = colData(int))
rm(ct, int)
sce <- sce[, which(!is.na(sce$infCellType))]
sce$infCellType <- droplevels(sce$infCellType)
### ### ###


### THEIR DATA - HYPEROXIA
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
samp <- 'P14_Hyperoxia'
sce <- hx[, which(hx$Sample == samp)]
### ### ###

# run CellChat
# prep SCE
assay(sce,'logcounts') <- CellChat::normalizeData(assay(sce,'counts'))
sce$samples <- factor(sce$Sample) # CellChat needs this name, specifically

cch <- createCellChat(object = sce, group.by = "CellType") # infCellType for our data
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB) # use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# set the used database in the object
cch@DB <- CellChatDB.use
cch <- subsetData(cch) # This step is necessary even if using the whole database (crazy)
cch <- identifyOverExpressedGenes(cch)
cch <- identifyOverExpressedInteractions(cch)
cch <- computeCommunProb(cch, type = "triMean")
cch <- filterCommunication(cch, min.cells = 10)
df.net <- subsetCommunication(cch)
cch <- computeCommunProbPathway(cch)
cch <- aggregateNet(cch)
# Compute the network centrality scores
cch <- netAnalysis_computeCentrality(cch, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

saveRDS(cch, file = 'data/cch_hx_P14Hyper.rds')

groupSize <- as.numeric(table(cch@idents))




