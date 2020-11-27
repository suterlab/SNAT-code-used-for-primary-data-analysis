library(Seurat)
library(tidyverse)
library(here)
library(cowplot)
source(here("Code", "utils.R"))

resDir <- here("Results", "SS2_P14_woFilter50k", "02-Analyses")
dir.create(resDir, recursive = TRUE)
figNumer <- 1

# Upgrade Seurat v2 object to v3 object
scData <- readRDS(here("Results", "SS2_P14_woFilter50k", "01-Seurat", "scData.rds"))
scData <- UpdateSeuratObject(scData)
scData <- RenameIdents(scData, "0"="tSC", "1"="mSC", "2"="pmSC", "3"="prol. SC")
scData <- RunUMAP(scData, dims = 1:7)
scData$Plate <- sub("_[[:alnum:]]+$", "", colnames(scData))
saveRDS(scData, file=file.path(resDir, "scData.rds"))

# Clustering plot
p1 <- DimPlot(scData, cols=cellCols, reduction = "umap", label=TRUE) +
  theme(legend.position = "none")
p2 <- DimPlot(scData, cols=cellCols, reduction = "tsne", label=TRUE) +
  theme(legend.position = "none")
p3 <- DimPlot(scData, reduction = "umap", label=TRUE, group.by="res.0.3") +
  theme(legend.position = "none")
p4 <- DimPlot(scData, reduction = "tsne", label=TRUE, group.by="res.0.3") +
  theme(legend.position = "none")
p <- plot_grid(p1, p2, p3, p4, labels="AUTO", ncol=2)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-ClusteringPlot.pdf")),
          p, ncol=2, nrow=2, base_asp = 1.2)
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
p <- DoHeatmap(scData, group.colors=cellCols[levels(Idents(scData))],
               features = top10$gene, raster=FALSE)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-pos_markers.pdf")),
          p, base_height = 6, base_width = 4)
figNumer <- figNumer + 1
