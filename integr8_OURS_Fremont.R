# integrate Jing's data (our data, 8 samples) with
# Fremont data (previously "ours"), PN14 Controls
require(Seurat)
cts <- Read10X('~/OneDrive - University of Southern California/lung_data/PN14_ctrl_filtered_feature_bc_matrix/')
so <- CreateSeuratObject(counts = cts)

# add genes
require(HDF5Array)
sce <- loadHDF5SummarizedExperiment('data/combined8filt_dbl/')
genes <- rownames(sce)
rm(sce)
toAdd <- genes[which(! genes %in% rownames(cts))]
if(length(toAdd) > 0){
  supp <- cts[seq_len(length(toAdd)), ]
  supp[,] <- 0
  rownames(supp) <- toAdd
  cts <- rbind(cts,supp)
  rm(supp)
}
cts <- cts[genes, ]

# make SCE
require(SingleCellExperiment)
fre <- SingleCellExperiment(assays = list(counts = cts))
fre$nCount <- so$nCount_RNA
fre$nFeature <- so$nFeature_RNA
rm(so, cts)

# null residuals
require(scry)
fre <- nullResiduals(fre, assay="counts", type="deviance")
fre$sample <- 'Fremont'

# add inferred cell type labels
fre_nat <- readRDS('data/integrated.rds')
cd <- colData(fre_nat)
cd <- cd[match(colnames(fre), cd$Cell_ID), ]
cd <- cd[,c('Cell_ID','infCellType','infCellType_conf')]
colData(fre) <- cbind(colData(fre), cd)
rm(fre_nat, cd)

# save as HDF5-backed SCE
saveHDF5SummarizedExperiment(fre, dir='data/tmp/fremont')



# combine with our data
# reset Fremont
rm(fre)
fre <- loadHDF5SummarizedExperiment('data/tmp/fremont/')
# load ours
ours <- loadHDF5SummarizedExperiment('data/combined8filt_dbl/')
# flatten doublet sub-tables
ours$dbl.clus.samp.class <- ours$dbl.clus.samp$class
ours$dbl.clus.samp.score <- ours$dbl.clus.samp$score
ours$dbl.noclus.samp.class <- ours$dbl.noclus.samp$class
ours$dbl.noclus.samp.score <- ours$dbl.noclus.samp$score
ours$dbl.clus.samp <- NULL
ours$dbl.noclus.samp <- NULL

# check:
stopifnot(all(rownames(fre)==rownames(ours)))
# make colData names the same (by adding NA columns)
for(n in names(colData(fre))){
  if(! n %in% names(colData(ours))){
    colData(ours)[[n]] <- NA
  }
}
for(n in names(colData(ours))){
  if(! n %in% names(colData(fre))){
    colData(fre)[[n]] <- NA
  }
}
colData(fre) <- colData(fre)[,names(colData(ours))]
# make new SCE with both datasets
sce <- SingleCellExperiment(
  assay = list(counts = cbind(assay(ours,'counts'), assay(fre,'counts')),
    binomial_deviance_residuals = cbind(assay(ours,'binomial_deviance_residuals'), assay(fre,'binomial_deviance_residuals'))),
  colData = rbind(colData(ours), colData(fre))
)
rm(ours,fre)

# sort by deviance
sce <- devianceFeatureSelection(sce, batch = factor(sce$sample), sorted = TRUE)


# fastMNN
require(batchelor)
sce <- fastMNN(sce, assay.type = 'binomial_deviance_residuals', batch = sce$orig.ident, subset.row = 1:2000) # top 2000 highest deviance genes


