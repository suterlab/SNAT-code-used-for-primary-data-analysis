library(Seurat)
library(tidyverse)
library(here)
library(cowplot)
library(stringr)
source(here("Code", "utils.R"))

resDir <- here("Results", "SS2_P60_woFilter50k", "02-Analyses")
dir.create(resDir, recursive = TRUE)
figNumer <- 1

# Upgrade Seurat v2 object to v3 object
scData <- readRDS(here("Results", "SS2_P60_woFilter50k", "01-Seurat", "scData.rds"))
scData <- UpdateSeuratObject(scData)
scData <- RunUMAP(scData, dims = 1:7)
scData$Plate <- str_extract(colnames(scData), "P60_(1|2)")
saveRDS(scData, file=file.path(resDir, "scData.rds"))

p <- DimPlot(scData, reduction = "tsne", label=TRUE) +
  theme(legend.position = "none")
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-ClusteringPlot_woLabels.pdf")),
          p, base_asp = 1.2)

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
p <- DoHeatmap(scData, features = top10$gene, raster=FALSE)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-pos_markers.pdf")),
          p, base_height = 6, base_width = 4)
figNumer <- figNumer + 1
