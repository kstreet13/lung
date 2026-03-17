# CellChat
# communication between fibroblasts and epithelial cells
# split by stages of development
# mesenchym ~= myofibroblast
# myofibroblast communication with AT1/AT2
# could compare our P14 to their P14, their P3/P7/P14

# differential signalling between control and hyperoxia

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


df.net <- subsetCommunication(cch)

cch <- computeCommunProbPathway(cch)

cch <- aggregateNet(cch)

groupSize <- as.numeric(table(cch@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cch@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cch@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cch, cluster.rows = TRUE, cluster.cols = TRUE, remove.isolate = TRUE)

netVisual_circle(cch@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",
                 idents.use = 'Myofibroblast')

netVisual_circle(cch@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",
                 idents.use = c('AT1','AT2 1','AT2 2'))


# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cch, sources.use = c('Myofibroblast'),
                 targets.use = c('AT1','AT2 1','AT2 2'),
                 remove.isolate = TRUE)

netVisual_bubble(cch, sources.use = c('AT1','AT2 1','AT2 2'),
                 targets.use = c('Myofibroblast'),
                 remove.isolate = TRUE)



# Compute the network centrality scores
cch <- netAnalysis_computeCentrality(cch, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cch, width = 12, height = 2.5, font.size = 10)

netAnalysis_signalingRole_scatter(cch)




