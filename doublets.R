# doublet detection
require(HDF5Array)
sce <- loadHDF5SummarizedExperiment('data/combined8filt_leiden')

require(scDblFinder)

# check with Jing: samples
# A vector [...] indicating to which sample each cell belongs.
# Here, a sample is understood as being processed independently. If omitted,
# doublets will be searched for with all cells together. If given, doublets will
# be searched for independently for each sample, which is preferable if they
# represent different captures.
# If your samples were multiplexed using cell hashes, what you want to give here
# are the different batches/wells (i.e. independent captures, since doublets
# cannot arise across them) rather than biological samples.

dbl.clus.samp <- scDblFinder(sce, clusters = 'leiden.r1', samples = 'sample',
                      includePCs = 25, returnType = 'scores', dbr.sd = 1)
dbl.noclus.samp <- scDblFinder(sce, samples = 'sample',
                       includePCs = 25, returnType = 'scores', dbr.sd = 1)

colData(sce)$dbl.clus.samp <- dbl.clus.samp
colData(sce)$dbl.noclus.samp <- dbl.noclus.samp

saveHDF5SummarizedExperiment(sce, dir = 'data/combined8filt_dbl')

# dbl.clus.nosamp <- scDblFinder(sce, clusters = 'leiden.r1',
#                              includePCs = 25, returnType = 'scores')
# dbl.noclus.nosamp <- scDblFinder(sce,
#                                includePCs = 25, returnType = 'scores')





ind <- sample(ncol(sce))
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = colorby(sce$dbl.noclus.samp$class[ind]), xlab = 'UMAP-1', ylab = 'UMAP-2')
legendby(sce$dbl.noclus.samp$class[ind])


