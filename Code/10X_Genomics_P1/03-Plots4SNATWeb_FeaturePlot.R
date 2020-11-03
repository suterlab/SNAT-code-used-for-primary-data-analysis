library(Seurat)
library(tidyverse)
library(here)
library(cowplot)
source(here("Code", "utils.R"))
im_size <- 6.945

resDir <- here("Results", "10x_P1", "03-Plots4SNATWeb_FeaturePlot")
dir.create(resDir, recursive = TRUE)

scData <- readRDS(here("Results", "10x_P1", "02-Analyses", "scData.rds"))

scData_all <- readRDS(here("Results", "10x_P1", "03-Plots4SNATWeb_SeuratAllGenes", "scData.rds"))
scData_all <- UpdateSeuratObject(scData_all)
scData_all <- scData_all[ ,colnames(scData)]
scData_all@reductions <- scData@reductions

g <- "Mbp"
for (g in rownames(scData_all)){
  message(g)
  p <- FeaturePlot(object = scData_all, features = g, reduction = "tsne",
                   cols = c("gray", "red"), pt.size=1) +
    panel_border(colour = "black", linetype = 1, remove = FALSE) +
    theme(legend.position="none", axis.title.x=element_blank(),
          axis.title.y=element_blank())

  save_plot(file.path(resDir, str_c(g, "_tSNE.pdf")),
            plot=p, base_height=im_size, base_width=im_size)
  save_plot(file.path(resDir, str_c(g, "_tSNE.png")),
            plot=p, base_height=im_size, base_width=im_size)
  save_plot(file.path(resDir, str_c(g, "_tSNE.eps")),
            plot=p, base_height=im_size, base_width=im_size)
}
