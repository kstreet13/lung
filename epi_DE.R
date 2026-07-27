################
# LOAD / RESET Epithelial SCE
################
require(HDF5Array)
require(SingleCellExperiment)
sce <- loadHDF5SummarizedExperiment('data/combined8reclus/')
sce <- sce[ ,which(sce$leiden.r1 %in% c(2,14,21,22,32,33,35))]
anno <- readRDS('data/epithelialANNO.rds')
sce$clus.epi <- anno$clus.epi
reducedDim(sce,'pca.epi') <- anno$pca.epi
reducedDim(sce,'umap2.epi') <- anno$umap2.epi
reducedDim(sce,'umap3.epi') <- anno$umap3.epi
rm(anno)


# DE between E4 and other AT2s
#  - what is E4?

##################
# MAKE SEURAT OBJECT
##################
require(Seurat); require(Matrix)
x <- sce
assayNames(x)[2] <- 'logcounts'
assay(x,'counts') <- realize(assay(x,'counts'))
assay(x,'counts') <- as(assay(x,'counts'), 'dgCMatrix')
assay(x,'logcounts') <- realize(assay(x,'logcounts'))
assay(x,'logcounts') <- as(assay(x,'logcounts'), 'dgCMatrix')
so <- as.Seurat(x)
rm(x)

###################
# DIFFERENTIAL EXPRESSION
###################
e4de <- FindMarkers(so, test.use = 'negbinom',
                    ident.1 = colnames(sce)[which(sce$clus.epi == 'E4')],
                    ident.2 = colnames(sce)[which(sce$clus.epi %in% c('E1','E2'))])
e4de$p_val_adj[e4de$p_val_adj == 0] <- min(e4de$p_val_adj[e4de$p_val_adj != 0]) / 2

# volcano
plot(e4de$avg_log2FC, -log10(e4de$p_val_adj))

tolabel <- which(e4de$avg_log2FC > 15 |
                   -log10(e4de$p_val_adj) > 70 |
                   e4de$avg_log2FC < -19 |
                   (-log10(e4de$p_val_adj) > 40 & e4de$avg_log2FC < -10)
)

plot(e4de$avg_log2FC, -log10(e4de$p_val_adj),
     col = rgb((e4de$p_val_adj < .01)/1.5, (abs(e4de$avg_log2FC) > 1)/1.5, 0,.25))
abline(v=0, lty=2)
text(e4de$avg_log2FC[tolabel], -log10(e4de$p_val_adj)[tolabel], rownames(e4de)[tolabel], cex=.75)

text(e4de$avg_log2FC, -log10(e4de$p_val_adj), rownames(e4de), cex=.5)








hornsde <- FindMarkers(so, test.use = 'negbinom',
                    ident.1 = colnames(sce)[which(sce$clus.epi == 'E5')],
                    ident.2 = colnames(sce)[which(sce$clus.epi == 'E6')])
hornsde$p_val_adj[hornsde$p_val_adj == 0] <- min(hornsde$p_val_adj[hornsde$p_val_adj != 0]) / 2


plot(hornsde$avg_log2FC, -log10(hornsde$p_val_adj),
     col = rgb((hornsde$p_val_adj < .01)/1.5, (abs(hornsde$avg_log2FC) > 1)/1.5, 0,.25))
abline(v=0, lty=2)
text(hornsde$avg_log2FC, -log10(hornsde$p_val_adj), rownames(hornsde), cex=.5)





at2de <- FindMarkers(so, test.use = 'negbinom',
                    ident.1 = colnames(sce)[which(sce$clus.epi == 'E3')],
                    ident.2 = colnames(sce)[which(sce$clus.epi %in% c('E1','E2'))])
at2de$p_val_adj[at2de$p_val_adj == 0] <- min(at2de$p_val_adj[at2de$p_val_adj != 0]) / 2

plot(at2de$avg_log2FC, -log10(at2de$p_val_adj),
     col = rgb((at2de$p_val_adj < .01)/1.5, (abs(at2de$avg_log2FC) > 1)/1.5, 0,.25))
abline(v=0, lty=2)
text(at2de$avg_log2FC, -log10(at2de$p_val_adj), rownames(at2de), cex=.5)

plot(-log10(at2de$p_val_adj), at2de$avg_log2FC,
     col = rgb((at2de$p_val_adj < .01)/1.5, (abs(at2de$avg_log2FC) > 1)/1.5, 0,.25))
abline(h=0, lty=2)
text(-log10(at2de$p_val_adj), at2de$avg_log2FC, rownames(at2de), cex=.5)







################
# GENE ONTOLOGY ANALYSIS
################
### functions
require(topGO)
require(limma)
addEntrezIDs <- function(df){
  require(org.Mm.eg.db)
  anno <- select(org.Mm.eg.db, keys=rownames(df), 
                 columns="ENTREZID", keytype="SYMBOL")
  stopifnot(all(rownames(df)==anno$SYMBOL))
  df$Entrez <- anno$ENTREZID
  df
}
topTerms <- function(res){
  stopifnot(any(c('Pathway','Term') %in% names(res)))
  type <- c('Pathway','Term')[which.max(c('Pathway','Term') %in% names(res))]
  ind <- which(res$p.adj < .05)
  if(length(ind)==0){
    ind <- 1:10
  }
  return(res[ind, ])
}
topTermsPlot <- function(res, topn = 10, ...){
  stopifnot(any(c('Pathway','Term') %in% names(res)))
  type <- c('Pathway','Term')[which.max(c('Pathway','Term') %in% names(res))]
  df <- res[topn:1, ]
  barplot(-log10(df$P.DE), horiz = TRUE, 
          xlab = '-log10 P-val', ylab = type, ...)
  t.x <- -log10(df$P.DE[nrow(df)])
  t.y <- seq(0.7, 0.7+1.2*(topn - 1), by = 1.2)
  text(t.x,t.y, df[[type]], pos = 2, cex=.8)
}



### Set comparison
de <- e4de
### e4de, hornsde, at2de

de <- addEntrezIDs(de)
univ <- select(org.Mm.eg.db, keys=rownames(sce), 
               columns="ENTREZID", keytype="SYMBOL")$ENTREZID

# *** split into up/down
deUP <- de[de$avg_log2FC > 0 & de$p_val_adj < .05, ]
deDN <- de[de$avg_log2FC < 0 & de$p_val_adj < .05, ]

goUP <- goana(deUP$Entrez, universe = univ, species = "Mm")
goUP <- goUP[order(goUP$P.DE), ]
goUP$p.adj <- p.adjust(goUP$P.DE, method = 'fdr')
goDN <- goana(deDN$Entrez, universe = univ, species = "Mm")
goDN <- goDN[order(goDN$P.DE), ]
goDN$p.adj <- p.adjust(goDN$P.DE, method = 'fdr')

topTerms(goUP)
topTerms(goDN)

topTermsPlot(goUP, main='Up in E5')
topTermsPlot(goDN, main='Down in E5')


# E4 is not noticeably lower quality
boxplot(sce$nCount_RNA ~ sce$clus.epi, log='y')
boxplot(sce$nFeature_RNA ~ sce$clus.epi, log='y')
