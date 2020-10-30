library(Seurat)
library(tidyverse)
library(here)
library(cowplot)
source(here("Code", "utils.R"))
# fixInNamespace("SingleExIPlot", pos="package:Seurat")
# jitter <- geom_jitter(height = 0, size = pt.size, alpha=0.5)

resDir <- here("Results", "10x_P60", "03-Plots4SNATWeb_ViolinPlot")
dir.create(resDir, recursive = TRUE)

scData <- readRDS(here("Results", "10x_P60", "02-Analyses", "scData.rds"))

scData_all <- readRDS(here("Results", "10x_P60", "03-Plots4SNATWeb_SeuratAllGenes", "scData.rds"))
scData_all <- UpdateSeuratObject(scData_all)
scData_all <- scData_all[ ,colnames(scData)]

scData@assays <- subset(scData_all, cells=colnames(scData))@assays

g <- "Mbp"
for (g in rownames(scData)){
  message(g)
  p <- VlnPlot(object = scData, features = g, cols = cellCols10xP60) +
    ylab("Expression (ln)") +
    theme(legend.position="none", axis.title.x=element_blank())

  save_plot(file.path(resDir, str_c(g, "_violins.pdf")),
            plot=p, base_height=4, base_width=6)
  save_plot(file.path(resDir, str_c(g, "_violins.png")),
            plot=p, base_height=4, base_width=6)
  save_plot(file.path(resDir, str_c(g, "_violins.eps")),
            plot=p, base_height=4, base_width=6)
}
