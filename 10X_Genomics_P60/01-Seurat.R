library(ezRun)
library(Seurat, lib.loc = "/home/gtan/R_Seurat2")
library(Matrix)
library(tidyverse)
library(here)
library(cowplot)

resDir <- here("Results", "10x_P60", "01-Seurat")
dir.create(resDir, recursive = TRUE)
figNumer <- 1

# Load the gene-cell count
cts <- list()
cts[[1]] <- Read10X(data.dir = "/srv/gstore/projects/p2497/CellRangerCount_24549_2018-05-23--00-43-09/P60_Sciatic_Nerve/outs/filtered_gene_bc_matrices/10X_Ref_Mouse_GRCm38.p5_20180305_Release_91/")
cts[[2]] <- Read10X(data.dir = "/srv/gstore/projects/p2497/CellRangerCount_26591_2018-05-11--14-23-14/Mouse_p60_1/outs/filtered_gene_bc_matrices/10X_Ref_Mouse_GRCm38.p5_20180305_Release_91/")
cts[[3]] <- Read10X(data.dir = "/srv/gstore/projects/p2497/CellRangerCount_26591_2018-05-11--14-23-14/Mouse_p60_2/outs/filtered_gene_bc_matrices/10X_Ref_Mouse_GRCm38.p5_20180305_Release_91/")
colnames(cts[[1]]) <- paste0("v1_", colnames(cts[[1]]))
colnames(cts[[2]]) <- paste0("v2_", colnames(cts[[2]]))
colnames(cts[[3]]) <- paste0("v3_", colnames(cts[[3]]))
cts <- cbind(cts[[1]], cts[[2]], cts[[3]])
scData <- CreateSeuratObject(raw.data = cts, min.cells = 5, project = "10X_P60")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = scData@data), value = TRUE)
percent.mito <- Matrix::colSums(scData@raw.data[mito.genes, ])/Matrix::colSums(scData@raw.data)
scData <- AddMetaData(object = scData, metadata = percent.mito,
                      col.name = "percent.mito")
p <- VlnPlot(object = scData, features.plot = c("nGene", "nUMI", "percent.mito"),
             nCol = 3)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-VlnPlot.pdf")),
          p, ncol=3)
figNumer <- figNumer + 1

scData <- FilterCells(object = scData, subset.names = c("nGene", "percent.mito"),
                      low.thresholds = c(500, -Inf), high.thresholds = c(5000, 0.25))
scData <- NormalizeData(object = scData, normalization.method = "LogNormalize",
                        scale.factor = 10000)

scData <- FindVariableGenes(object = scData, mean.function = ExpMean,
                            dispersion.function = LogVMR,
                            x.low.cutoff = 0.0125, x.high.cutoff = 3,
                            y.cutoff = 0.5)
length(x = scData@var.genes)
scData <- ScaleData(object = scData, vars.to.regress = c("nUMI", "percent.mito"))
scData <- RunPCA(object = scData, pc.genes = scData@var.genes, do.print = TRUE,
                 pcs.print = 1:5, genes.print = 5)
PrintPCA(object = scData, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
scData <- ProjectPCA(object = scData, do.print = FALSE)
scData <- FindClusters(object = scData, reduction.type = "pca", dims.use = 1:10,
                       resolution = 0.4, print.output = 0, prune.SNN = 0/100,
                       algorithm = 3, force.recalc=TRUE)
set.seed(100)
scData <- RunTSNE(object = scData, dims.use = 1:10, do.fast = TRUE, perplexity = 100)
TSNEPlot(object = scData, do.label = T)
saveRDS(scData, file=file.path(resDir, "scData.rds"))
