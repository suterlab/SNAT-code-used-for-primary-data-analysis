library(ezRun)
library(Seurat, lib.loc = "/home/gtan/R_Seurat2")
library(Matrix)
library(tidyverse)
library(here)
library(cowplot)

resDir <- here("Results", "SS2_P60_woFilter50k", "01-Seurat")
dir.create(resDir, recursive = TRUE)
figNumer <- 1

# Setup
min_genes=1500
max_genes=10000
project="p2497"
order="P60_merged"

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
rownames(cts) <- unname(gene_id2gene_name[rownames(cts)])

# Seurat workflow
scData <- CreateSeuratObject(raw.data = cts, min.cells = 0,
                             project = paste0(project, "_", order))
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
                      low.thresholds = c(min_genes, -Inf),
                      high.thresholds = c(max_genes, 0.25))
scData <- NormalizeData(object = scData, normalization.method = "LogNormalize",
                        scale.factor = 100000)

scData <- FindVariableGenes(object = scData)
length(x = scData@var.genes)
scData <- ScaleData(object = scData, vars.to.regress = c("nUMI","percent.mito"))
scData <- RunPCA(object = scData, pc.genes = scData@var.genes, do.print = TRUE,
                 pcs.print = 1:5, genes.print = 5)
scData <- ProjectPCA(object = scData, do.print = FALSE)
scData <- FindClusters(object = scData, reduction.type = "pca", dims.use = 1:7,
                       resolution = 0.3, print.output = 0, save.SNN = TRUE,
                       force.recalc=TRUE)
PrintFindClustersParams(object = scData)
set.seed(1)
scData <- RunTSNE(object = scData, dims.use = 1:7, do.fast = TRUE, perplexity = 30)
saveRDS(scData, file=file.path(resDir, "scData.rds"))

