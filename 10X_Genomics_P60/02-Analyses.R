library(Seurat)
library(tidyverse)
library(here)
library(cowplot)
source(here("Code", "utils.R"))

resDir <- here("Results", "10x_P60", "02-Analyses")
dir.create(resDir, recursive = TRUE)
figNumer <- 1

# Upgrade Seurat v2 object to v3 object
scData <- readRDS(here("Results", "10x_P60", "01-Seurat", "scData.rds"))
# scData <- readRDS("/home/gtan/analysis/p2497-DanielJorge/10X/10X_P60/scData.rds")
scData <- UpdateSeuratObject(scData)
p1 <- DimPlot(scData, reduction = "tsne", label=TRUE) +
  theme(legend.position = "none")
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-ClusteringPlotWOLabels.pdf")),
          p1, base_asp = 1.2)
figNumer <- figNumer + 1

## Fix some cells in cluster 7 to cluster 100
cellsInC7 <- WhichCells(scData, idents= 7)
tSNE_data <- as_tibble(Embeddings(scData, reduction="tsne"), rownames="cells")
tSNE_data <- tSNE_data %>% filter(tSNE_1 < 25 & tSNE_2 < -15 & cells %in% cellsInC7)
idents <- setNames(as.character(Idents(scData)), colnames(scData))
idents[tSNE_data$cells] <- "100"
Idents(scData) <- idents

scData <- RenameIdents(scData, "1"="EpC", "0"="EpC", "12"="EpC",
                       "10"="IC", "4"="IC",  "2"="PnC", "9"="PnC",
                       "5"="EnC", "11"="SC", "3"="EC1", "6"="EC2",
                       "7"="Per/VSMC", "8"="Per/VSMC", "100"="Per/EC*")
Idents(scData) <- factor(Idents(scData),
                         levels=c("SC", "EpC", "EnC", "PnC",
                                  "IC", "EC1", "EC2", "Per/EC*",
                                  "Per/VSMC"))
scData <- RunUMAP(scData, dims = 1:10)
scData$Plate <- sub("_[[:alnum:]]+$", "", colnames(scData)) %>%
  dplyr::recode("v1"="Run 1", "v2"="Run 2", "v3"="Run 3")
saveRDS(scData, file=file.path(resDir, "scData.rds"))

# Clustering plot
p1 <- DimPlot(scData, cols=cellCols10xP60, reduction = "umap", label=TRUE) +
  theme(legend.position = "none")
p2 <- DimPlot(scData, cols=cellCols10xP60, reduction = "tsne", label=TRUE) +
  theme(legend.position = "none")
p <- plot_grid(p1, p2, labels="AUTO", ncol=2)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-ClusteringPlot.pdf")),
          p, ncol=2, base_asp = 1.2)
figNumer <- figNumer + 1

# QC plots
p1 <- FeaturePlot(object=scData, features="nFeature_RNA",
                  reduction="umap", cols=c("yellow", "red"))
p2 <- FeaturePlot(object=scData, features="nCount_RNA",
                  reduction="umap", cols=c("yellow", "red"))
p3 <- FeaturePlot(object=scData, features="percent.mito",
                  reduction="umap", cols=c("yellow", "red"))
p <- plot_grid(p1, p2, p3, ncol=3)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-QCPlot_umap.pdf")),
          p, ncol=3)
p1 <- FeaturePlot(object=scData, features="nFeature_RNA",
                  reduction="tsne", cols=c("yellow", "red"))
p2 <- FeaturePlot(object=scData, features="nCount_RNA",
                  reduction="tsne", cols=c("yellow", "red"))
p3 <- FeaturePlot(object=scData, features="percent.mito",
                  reduction="tsne", cols=c("yellow", "red"))
p <- plot_grid(p1, p2, p3, ncol=3)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-QCPlot_tsne.pdf")),
          p, ncol=3)
figNumer <- figNumer + 1

# Plate origin
p1 <- DimPlot(scData, reduction = "umap", group.by="Plate")
p2 <- DimPlot(scData, reduction = "tsne", group.by="Plate")
p <- plot_grid(p1, p2, labels="AUTO", ncol=2)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-PlateOrigin.pdf")),
          p, ncol=2, base_asp = 1.4)
figNumer <- figNumer + 1

# markers, heatmap
markers <- FindAllMarkers(scData, only.pos=TRUE, return.thresh=0.01)
saveRDS(markers, file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                          "-pos_markers.rds")))
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
p <- DoHeatmap(scData, group.colors=cellCols10xP60[levels(Idents(scData))],
               features = top10$gene, raster=FALSE)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-pos_markers.pdf")),
          p, base_height = 14, base_asp = 0.6)
figNumer <- figNumer + 1
