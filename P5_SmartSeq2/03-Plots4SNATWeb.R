library(Seurat)
library(tidyverse)
library(here)
library(cowplot)
source(here("Code", "utils.R"))

resDir <- here("Results", "SS2_P5_woFilter50k", "03-Plots4SNATWeb")
dir.create(resDir, recursive = TRUE)

im_size <- 6.945
scData <- readRDS(here("Results", "SS2_P5_woFilter50k", "02-Analyses", "scData.rds"))

g <- rownames(scData)[1]
for (g in rownames(scData)){
  message(g)
  p <- FeaturePlot(object = scData, features = g, reduction = "tsne",
                   cols = c("gray", "red"), pt.size = 2) +
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

