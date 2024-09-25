
# compare their data with "our" data
# PN14 Ctrl (ours) vs. hyperoxia (theirs, sce)


require(SingleCellExperiment)
hyperoxia <- readRDS('~/OneDrive - University of Southern California/lung_data/sce.rds')

require(DropletUtils)
ctrl <- read10xCounts('~/OneDrive - University of Southern California/lung_data/PN14_Ctrl_filtered_feature_bc_matrix/', sample.names = "PN14ctrl")

