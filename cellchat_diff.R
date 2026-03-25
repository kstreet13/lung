# differential signalling between control and hyperoxia

norm <- readRDS('data/cch_hx_P14Norm.rds')
hyp <- readRDS('data/cch_hx_P14Hyper.rds')


object.list <- list(norm = norm, hyp = hyp)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# number/strength of interactions
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# differential interactions
# not helpful, too busy
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# RED = higher in Hyperoxia
# BLUE = higher in Normoxia
netVisual_heatmap(cellchat, measure = "weight")

# incoming/outgoing scatter 
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg1 <- netAnalysis_signalingRole_scatter(object.list[[1]], title = names(object.list)[1], weight.MinMax = weight.MinMax)
gg2 <- netAnalysis_signalingRole_scatter(object.list[[2]], title = names(object.list)[2], weight.MinMax = weight.MinMax)
gg1 + gg2

# identify the specific signaling changes for given cell type
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Myofibroblast")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "AT1")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "AT2 1")
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "AT2 2")
patchwork::wrap_plots(plots = list(gg1,gg2,gg3,gg4), ncol=2)



cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

# cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural")
# netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)




gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "hyp"
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "hyp",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "norm",ligand.logFC = -0.05, receptor.logFC = NULL)

# up in hyperoxia
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)

# up in normoxia
computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)

