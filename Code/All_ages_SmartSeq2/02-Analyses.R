library(Seurat)
library(tidyverse)
library(here)
library(cowplot)
library(stringr)
source(here("Code", "utils.R"))

resDir <- here("Results", "SS2_P1P5P14P60_merged", "02-Analyses")
dir.create(resDir, recursive = TRUE)
figNumer <- 1

# Upgrade Seurat v2 object to v3 object
scData <- readRDS(here("Results", "SS2_P1P5P14P60_merged", "01-Seurat", "scData.rds"))
#scData <- readRDS("/srv/GT/analysis/p2497/seurat_P1P5P14P60_merged_scData_vplots.rds")
scData <- UpdateSeuratObject(scData)
scData <- RenameIdents(scData, "0"="iSC", "1"="nm(R)SC", "2"="mSC cluster 2", "3"="tSC",
                       "4"="mSC cluster 1", "5"="pmSC", "6"="prol. SC", "7"="mSC cluster 3")
Idents(scData) <- factor(Idents(scData),
                         levels=c("prol. SC", "iSC", "pmSC", "mSC cluster 1", "mSC cluster 2",
                                  "mSC cluster 3", "tSC", "nm(R)SC"))
scData <- RunUMAP(scData, dims = 1:20)
scData$Plate <- str_extract(colnames(scData), "P\\d+_(1|2)")
scData$timepoint <- sub("_(1|2)$", "", scData$Plate)
saveRDS(scData, file=file.path(resDir, "scData.rds"))

# Clustering plot
p1 <- DimPlot(scData, cols=cellCols, reduction = "umap", label=TRUE) +
  theme(legend.position = "none")
p2 <- DimPlot(scData, cols=cellCols, reduction = "tsne", label=TRUE) +
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
                  reduction="umap", cols=c("yellow", "red")) +
  ggtitle("fraction mito")
p <- plot_grid(p1, p2, p3, ncol=3)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-QCPlot_umap.pdf")),
          p, ncol=3)
p1 <- FeaturePlot(object=scData, features="nFeature_RNA",
                  reduction="tsne", cols=c("yellow", "red"))
p2 <- FeaturePlot(object=scData, features="nCount_RNA",
                  reduction="tsne", cols=c("yellow", "red"))
p3 <- FeaturePlot(object=scData, features="percent.mito",
                  reduction="tsne", cols=c("yellow", "red")) +
  ggtitle("fraction mito")
p <- plot_grid(p1, p2, p3, ncol=3)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-QCPlot_tsne.pdf")),
          p, ncol=3)
figNumer <- figNumer + 1

# Plate origin
p1 <- DimPlot(scData, cols="Paired", reduction = "umap", group.by="Plate")
p2 <- DimPlot(scData, cols="Paired", reduction = "tsne", group.by="Plate")
p <- plot_grid(p1, p2, labels="AUTO", ncol=2)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-PlateOrigin.pdf")),
          p, ncol=2, base_asp = 1.4)

# timecourse
p1 <- DimPlot(scData, reduction = "umap", group.by="timepoint")
p2 <- DimPlot(scData, reduction = "tsne", group.by="timepoint")
p <- plot_grid(p1, p2, labels="AUTO", ncol=2)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-timepoint.pdf")),
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
          p, base_height = 10, base_width = 6)
figNumer <- figNumer + 1
