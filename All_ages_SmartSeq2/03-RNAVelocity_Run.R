# Run RNA velocity on SS2 merged samples
library(reticulate)
library(ezRun)
library(here)
use_condaenv("python3", conda = "/usr/local/ngseq/miniconda3/bin/conda",
             required = FALSE)
py_discover_config()
cwd <- getwd()
resDir <- here("Results", "SS2_P1P5P14P60_merged", "03-RNAVelocity_Run")
dir.create(resDir, recursive = TRUE)

## Prepare bam files
bamFns <- c("P1_1"="projects/p2497/SCCountsApp_27067_2018-07-02--18-01-13/P1_1.bam",
            "P1_2"="projects/p2497/SCCountsApp_27067_2018-07-02--18-01-13/P1_2.bam",
            "P5_1"="projects/p2497/SCCountsApp_27068_2018-06-10--00-28-19/P5_1.bam",
            "P5_2"="projects/p2497/SCCountsApp_27068_2018-06-10--00-28-19/P5_2.bam",
            "P14_1"="projects/p2497/SCCountsApp_27068_2018-06-10--00-28-19/P14_1.bam",
            "P14_2"="projects/p2497/SCCountsApp_27068_2018-06-10--00-28-19/P14_2.bam",
            "P60_1"="projects/p2497/SCCountsApp_24885_2018-07-03--01-58-26/20180221.A-SiCSeq_SCs_P5_E1_unmapped.bam",
            "P60_2"="projects/p2497/SCCountsApp_26718_2018-05-24--20-38-34/P60_2.bam"
            )
bamFns <- setNames(file.path("/srv/gstore", bamFns), names(bamFns))

### use local scratch; nfs will be too slow
setwdNew("/scratch/gtan/p2497-DanielJorge/velocity/SS2_merged/Bams")
bams <- lapply(bamFns, splitBamByRG, mc.cores=8L)
bams <- unlist(bams)

## Run velocyto
cmd <- paste("velocyto run-smartseq2", "-o SS2_merged_velocity",
             "-e SS2_merged", "-v",
             paste(bams, collapse=" "),
             "/srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26/Genes/genes.gtf")
ezSystem(cmd)

file.copy("SS2_merged_velocity", resDir, recursive = TRUE)
