library(Seurat)
library(tidyverse)
library(here)
library(cowplot)
library(fastICA)
library(mgcv)
source(here("Code", "utils.R"))
im_size <- 5

resDir <- here("Results", "SS2_P1P5P14P60_merged", "06-Plots4SNATWeb_Trajectory")
dir.create(resDir, recursive = TRUE)

scData <- readRDS(here("Results", "SS2_P1P5P14P60_merged", "02-Analyses", "scData.rds"))
scData_all <- readRDS(here("Results", "SS2_P1P5P14P60_merged", "06-Plots4SNATWeb_SeuratAllGenes", "scData.rds"))
scData_all <- UpdateSeuratObject(scData_all)

### smoothing spline: myelinating trajectory
scData_m <- subset(scData, idents = c("prol. SC", "iSC", "pmSC",
                                      "mSC cluster 1", "mSC cluster 2",
                                      "mSC cluster 3"))
set.seed(1)
log_cd <- GetAssayData(scData_m)
a <- fastICA(log_cd, 2, alg.typ = "deflation", fun = "logcosh", alpha = 1,
             method = "C", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)
dt <- data.frame(row.names=colnames(scData_m),
                 clusters=Idents(scData_m),
                 pseudotime=a$A[1, ])
pdata <- data.frame(pseudotime=-200:200/200)

scData_m@assays <- subset(scData_all, cells=colnames(scData_m))@assays
g <- "Mbp"
for (g in rownames(scData_m)){
  message(g)
  toPlot <- bind_cols(dt, gene_exp=GetAssayData(scData_m)[g, ])
  gam_model <- gam(gene_exp~s(pseudotime, k=5), data=toPlot)
  pv <- predict(gam_model, pdata, type="response", se=TRUE)
  toPlot2 <- bind_cols(pdata, fit=pv$fit, fit_high=pv$fit+1.96*pv$se.fit,
                       fit_low=pv$fit-1.96*pv$se.fit)

  p <- ggplot(toPlot, aes(-pseudotime, gene_exp)) +
    geom_point(aes(colour=clusters)) +
    scale_color_manual(values=cellCols_SS2Merged) +
    xlab("myelinating trajectory pseudotime") +
    ylab("Expression (ln)") +
    ggtitle(g) +
    geom_line(aes(-pseudotime, fit), data=toPlot2) +
    geom_line(aes(-pseudotime, fit_high), data=toPlot2, linetype = 2) +
    geom_line(aes(-pseudotime, fit_low), data=toPlot2, linetype = 2) +
    theme_cowplot() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    xlim(range(-toPlot$pseudotime)) + ylim(range(toPlot$gene_exp))

  save_plot(filename=file.path(resDir, str_c(g, "_trajectory_myelinating.pdf")),
            p, base_height = im_size, base_width = im_size)
  save_plot(filename=file.path(resDir, str_c(g, "_trajectory_myelinating.eps")),
            p, base_height = im_size, base_width = im_size)
  save_plot(filename=file.path(resDir, str_c(g, "_trajectory_myelinating.png")),
            p, base_height = im_size, base_width = im_size)
}

### smoothing spline: non-myelinating trajectory
scData_m <- subset(scData, idents = c("prol. SC", "iSC", "tSC", "nm(R)SC"))
set.seed(1)
log_cd <- GetAssayData(scData_m)
a <- fastICA(log_cd, 2, alg.typ = "deflation", fun = "logcosh", alpha = 1,
             method = "C", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)
dt <- data.frame(row.names=colnames(scData_m),
                 clusters=Idents(scData_m),
                 pseudotime=a$A[1, ])
pdata <- data.frame(pseudotime=-200:200/200)

scData_m@assays <- subset(scData_all, cells=colnames(scData_m))@assays
g <- "Mbp"
for (g in rownames(scData_m)){
  message(g)
  toPlot <- bind_cols(dt, gene_exp=GetAssayData(scData_m)[g, ])
  gam_model <- gam(gene_exp~s(pseudotime, k=5), data=toPlot)
  pv <- predict(gam_model, pdata, type="response", se=TRUE)
  toPlot2 <- bind_cols(pdata, fit=pv$fit, fit_high=pv$fit+1.96*pv$se.fit,
                       fit_low=pv$fit-1.96*pv$se.fit)

  p <- ggplot(toPlot, aes(-pseudotime, gene_exp)) +
    geom_point(aes(colour=clusters)) +
    scale_color_manual(values=cellCols_SS2Merged) +
    xlab("non-myelinating trajectory pseudotime") +
    ylab("Expression (ln)") +
    ggtitle(g) +
    geom_line(aes(-pseudotime, fit), data=toPlot2) +
    geom_line(aes(-pseudotime, fit_high), data=toPlot2, linetype = 2) +
    geom_line(aes(-pseudotime, fit_low), data=toPlot2, linetype = 2) +
    theme_cowplot() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    xlim(range(-toPlot$pseudotime)) + ylim(range(toPlot$gene_exp))

  save_plot(filename=file.path(resDir, str_c(g, "_trajectory_non-myelinating.pdf")),
            p, base_height = im_size, base_width = im_size)
  save_plot(filename=file.path(resDir, str_c(g, "_trajectory_non-myelinating.eps")),
            p, base_height = im_size, base_width = im_size)
  save_plot(filename=file.path(resDir, str_c(g, "_trajectory_non-myelinating.png")),
            p, base_height = im_size, base_width = im_size)
}
