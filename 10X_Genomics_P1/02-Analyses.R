library(Seurat)
library(tidyverse)
library(here)
library(cowplot)
source(here("Code", "utils.R"))

resDir <- here("Results", "10x_P1", "02-Analyses")
dir.create(resDir, recursive = TRUE)
figNumer <- 1

# Upgrade Seurat v2 object to v3 object
scData <- readRDS(here("Results", "10x_P1", "01-Seurat", "scData.rds"))
scData <- readRDS("/srv/GT/analysis/p2497/10X_P1_merged_scData_vplots.rds")
scData <- UpdateSeuratObject(scData)
scData <- RenameIdents(scData, "1"="pmSC", "4"="iSC", "9"="prol. SC",
                       "11"="IC", "5"="IC", "6"="EC", "7"="prol. Fb",
                       "2"="EnC", "10"="FbRel*", "0"="EpC", "3"="PnC",
                       "8"="Per/VSMC")
Idents(scData) <- factor(Idents(scData),
                         levels=c("prol. SC", "iSC", "pmSC", "prol. Fb",
                                  "EpC", "EnC", "PnC", "FbRel*",
                                  "Per/VSMC", "EC", "IC"))
scData <- RunUMAP(scData, dims = 1:10)
scData$Plate <- sub("_[[:alnum:]]+$", "", colnames(scData)) %>%
  dplyr::recode("v1"="Run 1", "v2"="Run 2", "v3"="Run 3")
saveRDS(scData, file=file.path(resDir, "scData.rds"))

# Clustering plot
p1 <- DimPlot(scData, cols=cellCols10xP1, reduction = "umap", label=TRUE) +
  theme(legend.position = "none")
p2 <- DimPlot(scData, cols=cellCols10xP1, reduction = "tsne", label=TRUE) +
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
p <- DoHeatmap(scData, group.colors=cellCols10xP1[levels(Idents(scData))],
               features = top10$gene, raster=FALSE)
save_plot(filename=file.path(resDir, paste0(formatC(figNumer, width=2, flag="0"),
                                            "-pos_markers.pdf")),
          p, base_height = 14, base_asp = 0.6)
figNumer <- figNumer + 1
