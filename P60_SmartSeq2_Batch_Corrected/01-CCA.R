library(ezRun)
library(Seurat, lib.loc = "/home/gtan/R_Seurat2")
library(Matrix)
library(tidyverse)
library(here)
library(cowplot)

resDir <- here("Results", "SS2_P60_BatchCorrection_woFilter50k", "01-CCA")
dir.create(resDir, recursive = TRUE)
figNumer <- 1

# Setup
min_genes <- 1500
max_genes <- 10000
project <- "p2497"
order <- "P60_merged"

# Load the gene-cell count
loc <- "SCCountsApp_24885_2018-07-03--01-58-26"
cts1 <- readSCMM(file.path("/srv/gstore/projects/p2497", loc, "20180221.A-SiCSeq_SCs_P5_E1_unmapped-counts.mtx"))
colnames(cts1) <- paste0("P60_1_", colnames(cts1))
loc <- "SCCountsApp_26718_2018-05-24--20-38-34"
cts2 <- readSCMM(file.path("/srv/gstore/projects/p2497", loc, "P60_2-counts.mtx"))
stopifnot(identical(rownames(cts1), rownames(cts2)))
cts <- cbind(cts1, cts2)

refBuild <- read_tsv(file.path("/srv/gstore/projects/p2497/", loc,
                               "/parameters.tsv")) %>%
  filter(sushi_app == "refBuild") %>% pull("SCCountsApp")
seqAnno <- read_tsv(file.path("/srv/GT/reference", refBuild,
                              "Genes/genes_annotation_byGene.txt"))
gene_id2gene_name <- setNames(scater::uniquifyFeatureNames(seqAnno$gene_id, seqAnno$gene_name),
                              seqAnno$gene_id)
rownames(cts1) <- unname(gene_id2gene_name[rownames(cts1)])
rownames(cts2) <- unname(gene_id2gene_name[rownames(cts2)])

# CCA correction
scData1 <- CreateSeuratObject(raw.data = cts1, project = "P60_1", min.cells = 0)
mito.genes <- grep(pattern = "^mt-", x = rownames(x = scData1@data), value = TRUE)
percent.mito <- Matrix::colSums(scData1@raw.data[mito.genes, ])/Matrix::colSums(scData1@raw.data)
scData1 <- AddMetaData(object = scData1, metadata = percent.mito, col.name = "percent.mito")
scData1@meta.data$stim <- "P60_1"
scData1 <- FilterCells(scData1, subset.names = c("nGene", "percent.mito"),
                       low.thresholds = c(min_genes, -Inf),
                       high.thresholds = c(max_genes, 0.25))
scData1 <- NormalizeData(scData1)
scData1 <- ScaleData(scData1, display.progress = F)

scData2 <- CreateSeuratObject(raw.data = cts2, project = "P60_2", min.cells = 0)
mito.genes <- grep(pattern = "^mt-", x = rownames(x = scData2@data), value = TRUE)
percent.mito <- Matrix::colSums(scData2@raw.data[mito.genes, ])/Matrix::colSums(scData2@raw.data)
scData2 <- AddMetaData(object = scData2, metadata = percent.mito, col.name = "percent.mito")
scData2@meta.data$stim <- "P60_2"
scData2 <- FilterCells(scData2, subset.names = c("nGene", "percent.mito"),
                       low.thresholds = c(min_genes, -Inf),
                       high.thresholds = c(max_genes, 0.25))

scData2 <- NormalizeData(scData2)
scData2 <- ScaleData(scData2, display.progress = F)

scData1 <- FindVariableGenes(scData1, do.plot = F)
scData2 <- FindVariableGenes(scData2, do.plot = F)
g.1 <- head(rownames(scData1@hvg.info), 1000)
g.2 <- head(rownames(scData2@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(scData1@scale.data))
genes.use <- intersect(genes.use, rownames(scData2@scale.data))

scData <- RunCCA(scData1, scData2, genes.use = genes.use, num.cc = 30)
p1 <- DimPlot(object = scData, reduction.use = "cca", group.by = "stim",
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = scData, features.plot = "CC1", group.by = "stim",
              do.return = TRUE)
plot_grid(p1, p2)
p3 <- MetageneBicorPlot(scData, grouping.var = "stim", dims.eval = 1:30,
                        display.progress = FALSE)
scData <- AlignSubspace(scData, reduction.type = "cca", grouping.var = "stim",
                        dims.align = 1:15)

p1 <- VlnPlot(object = scData, features.plot = "ACC1", group.by = "stim",
              do.return = TRUE)
p2 <- VlnPlot(object = scData, features.plot = "ACC2", group.by = "stim",
              do.return = TRUE)
plot_grid(p1, p2)

scData <- RunTSNE(scData, reduction.use = "cca.aligned", dims.use = 1:12,
                  do.fast = T)
scData <- FindClusters(scData, reduction.type = "cca.aligned",
                       resolution = 0.5, dims.use = 1:12)
p1 <- TSNEPlot(scData, do.return = T, group.by = "stim")
p2 <- TSNEPlot(scData, do.label = T, do.return = T)
plot_grid(p1, p2)

saveRDS(scData, file=file.path(resDir, "scData.rds"))
