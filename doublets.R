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

dblScr <- scDblFinder(sce, clusters = 'leiden.r1', samples = 'sample',
                      includePCs = 25, returnType = 'scores')
dblScr2 <- scDblFinder(sce, samples = 'sample',
                       includePCs = 25, returnType = 'scores')


ind <- sample(ncol(sce))
plot(reducedDim(sce,'umap')[ind,], asp=1, cex=.25, col = colorby(dblScr2$class[ind]), xlab = 'UMAP-1', ylab = 'UMAP-2')
legendby(dblScr2$class)


