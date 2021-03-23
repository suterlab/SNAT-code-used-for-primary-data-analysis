library(here)
library(tidyverse)
library(velocyto.R)
library(Seurat)
library(cowplot)
source(here("Code", "utils.R"))

resDir <- here("Results", "SS2_P1P5P14P60_merged", "03-RNAVelocity_Analyses")
dir.create(resDir, recursive = TRUE)

# Data loading
ldat <- read.loom.matrices(here("Results", "SS2_P1P5P14P60_merged", "03-RNAVelocity_Run",
                                "SS2_merged_velocity", "SS2_merged.loom"))
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub("\\.bam$", "", gsub(".*:","",colnames(x)))
  colnames(x) <- sub("^N", "P60_1_N", colnames(x))
  x
})

scData <- readRDS(here("Results", "SS2_P1P5P14P60_merged", "02-Analyses", "scData.rds"))
emb <- Embeddings(object = scData, reduction = "tsne")
cell.colors <- setNames(cellCols_SS2Merged[as.character(Idents(scData))], colnames(scData))
plot(emb[, "tSNE_1"], emb[, "tSNE_2"], col=cell.colors, pch=16)

emat_cutoff <- 5
nmat_cutoff <- 1
smat_cutoff <- 0.5
fit.quantile <- 0.05
kCells <- 5
cell.cex <- 1
arrow.scale <- 3
show.grid.flow <- FALSE
grid.n <- 20
cell.alpha <- 0.7

# Gene filtering
hist(log10(rowSums(ldat$spliced[intersect(rownames(ldat$spliced),
                                          rownames(scData)), ])+1), col='wheat',
     xlab='log10[ number of reads + 1]',
     main='number of reads per gene')

# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning

# restrict to cells that from filtering upstream
emat <- emat[ , rownames(emb)]
nmat <- nmat[ , rownames(emb)]
smat <- smat[ , rownames(emb)]

# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat, cell.colors,
                                           min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat, cell.colors,
                                           min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat, cell.colors,
                                           min.max.cluster.average = 0.5)
# look at the resulting gene set
length(intersect(rownames(emat),rownames(nmat)))
# and if we use spanning reads (smat)
length(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))

### Velocity estimates variant 1
rvel.qf <- gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells = 5,
                                            fit.quantile = 0.05)
CNEr:::savefig(file.path(resDir, "Variant1_PCA"), type="pdf", colormode="srgb",
               height=20, width=20)
pca.velocity.plot(rvel.qf, nPcs=5, plot.cols=2,
                  cell.colors=ac(cell.colors,alpha=0.7),
                  cex=1.2, pcount=0.1, pc.multipliers=c(1,-1,-1,-1,-1),
                  show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                  grid.n=20, arrow.lwd=2)
dev.off()

CNEr:::savefig(file.path(resDir, "Variant1_tSNE"), type="pdf", colormode="srgb",
               height=15, width=17)
show.velocity.on.embedding.cor(emb, rvel.qf, n=100, scale='sqrt',
                               cell.colors=ac(cell.colors,
                                              alpha=cell.alpha),
                               cex=cell.cex, arrow.scale=arrow.scale,
                               arrow.lwd=1)
legend("topleft", legend=names(cellCols_SS2Merged),
       pch=19, col=cellCols_SS2Merged)

show.velocity.on.embedding.cor(emb, rvel.qf, n=100, scale='sqrt',
                               cell.colors=ac(cell.colors,
                                              alpha=cell.alpha),
                               cex=cell.cex, arrow.scale=arrow.scale,
                               show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=grid.n,
                               arrow.lwd=2)
legend("topleft", legend=names(cellCols_SS2Merged),
       pch=19, col=cellCols_SS2Merged)
dev.off()

### Velocity estimates variant 2
rvel <- gene.relative.velocity.estimates(emat,nmat,smat=smat,
                                         kCells=kCells,
                                         fit.quantile=fit.quantile,
                                         diagonal.quantiles=TRUE, n.cores=8L)
CNEr:::savefig(file.path(resDir, "Variant2_PCA"), type="pdf", colormode="srgb",
               height=20, width=20)
pca.velocity.plot(rvel, nPcs=5, plot.cols=2,
                  cell.colors=ac(cell.colors,alpha=0.7),
                  cex=1.2, pcount=0.1, pc.multipliers=c(1,-1,1,1,1),
                  show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                  grid.n=grid.n, arrow.lwd=2, n.cores=8)
dev.off()

CNEr:::savefig(file.path(resDir, "Variant2_tSNE"), type="pdf", colormode="srgb",
               height=15, width=17)
show.velocity.on.embedding.cor(emb, rvel, n=100, scale='sqrt',
                               cell.colors=ac(cell.colors,
                                              alpha=cell.alpha),
                               cex=cell.cex, arrow.scale=arrow.scale,
                               arrow.lwd=1, n.cores=8L)
legend("topleft", legend=names(cellCols_SS2Merged),
       pch=19, col=cellCols_SS2Merged)
show.velocity.on.embedding.cor(emb, rvel, n=100, scale='sqrt',
                               cell.colors=ac(cell.colors,
                                              alpha=0.7),
                               cex=cell.cex, arrow.scale=arrow.scale,
                               show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                               grid.n=grid.n, arrow.lwd=2, n.cores=8L)
legend("topleft", legend=names(cellCols_SS2Merged),
       pch=19, col=cellCols_SS2Merged)
dev.off()
