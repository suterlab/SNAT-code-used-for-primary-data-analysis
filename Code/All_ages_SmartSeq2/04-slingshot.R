library(here)
library(Seurat)
library(slingshot)
library(cowplot)
library(tidyverse)
source(here("Code", "utils.R"))

resDir <- here("Results", "SS2_P1P5P14P60_merged", "04-slingshot")
dir.create(resDir, recursive = TRUE)

scData <- readRDS(here("Results", "SS2_P1P5P14P60_merged", "02-Analyses", "scData.rds"))
sce <- as.SingleCellExperiment(scData)
sce <- slingshot(sce, clusterLabels = 'ident', reducedDim = 'TSNE',
                 start.clus = "prol. SC")
saveRDS(sce, file.path(resDir, "sce.rds"))

sce <- readRDS(file.path(resDir, "sce.rds"))
CNEr:::savefig(filename = file.path(resDir, "slingshot"), type="pdf",
               colormodel = "srgb", height = 6*2.54, width=12*2.54)
par(mfrow=c(1,2))
plot(reducedDims(sce)$TSNE, col = cellCols[as.character(sce$ident)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=3, type = 'lineages')
plot(reducedDims(sce)$TSNE, col = cellCols[as.character(sce$ident)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=3)
dev.off()

