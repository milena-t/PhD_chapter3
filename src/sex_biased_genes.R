# install DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2") # compiled from source 
# install tidyverse (ggplot2 is in tidyverse)
install.packages("tidyverse")
library("DESeq2")
library(tidyverse)
library(ggplot2)

## load data
setwd("/Users/miltr339/work/chapter3/DE_analysis")
Cmac_path <- "Cmac_gene_counts_short_headers.txt"
Csep_path <- "Csep_gene_counts_short_headers.txt"
Tcas_path <- "Tcas_gene_counts_short_headers.txt"

## metadata relative paths
Cmac_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Cmac_full_SRR_list.csv'
Csep_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Csep_full_SRR_list.csv'
Tcas_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Tcas_full_SRR_list.csv'

## read data for DE
count_data <- read.csv(Cmac_path, header = TRUE, sep = "\t",comment.char = "#")
count_data <- count_data %>% select(-one_of('Chr','Start','End','Strand','Length')) 
DE_metadata_all <- read.csv(Cmac_metadata, header = TRUE, sep = ",",comment.char = "#")
# only include metadata when the samples are in the counts
meta_data <- DE_metadata_all[which((DE_metadata_all$Run %in% names(DE_counts_raw))==TRUE), ]

## set up DESeq data structures
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = meta_data,
                              design = organ ~ sex)
