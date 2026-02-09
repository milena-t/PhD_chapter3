### Everything is mostly based off of this blog post:
# https://ashleyschwartz.com/posts/2023/05/deseq2-tutorial
# and i use vst for count normalization: https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/varianceStabilizingTransformation


# install DESeq2
if (FALSE){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2") # compiled from source 
  # install tidyverse (ggplot2 is in tidyverse)
  install.packages("tidyverse") 
  install.packages("gtable") 
}
library("DESeq2")
library(tidyverse)
library(ggplot2)

#### load data
setwd("/Users/miltr339/work/chapter3/DE_analysis")
Cmac_path <- "Cmac_gene_counts_short_headers.txt"
Csep_path <- "Csep_gene_counts_short_headers.txt"
Tcas_path <- "Tcas_gene_counts_short_headers.txt"
## normalized counts outfiles
Cmac_path_out <- "Cmac_gene_counts_short_headers_vst.txt"
Csep_path_out <- "Csep_gene_counts_short_headers_vst.txt"
Tcas_path_out <- "Tcas_gene_counts_short_headers_vst.txt"
## log-fold-change outfiles
Cmac_path_out_lfc <- "Cmac_gene_counts_short_headers_lfc.txt"
Csep_path_out_lfc <- "Csep_gene_counts_short_headers_lfc.txt"
Tcas_path_out_lfc <- "Tcas_gene_counts_short_headers_lfc.txt"

#### metadata relative paths
Cmac_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Cmac_full_SRR_list.csv'
Csep_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Csep_full_SRR_list.csv'
Tcas_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Tcas_full_SRR_list.csv'

#### do one species of the three
if (FALSE){
  data_path <- Cmac_path
  out_path_vst <- Cmac_path_out
  out_path_lfc <- Cmac_path_out_lfc
  metadata <- Cmac_metadata
}
if (TRUE){
  data_path <- Csep_path
  out_path_vst <- Csep_path_out
  out_path_lfc <- Csep_path_out_lfc
  metadata <- Csep_metadata
}
if (FALSE){
  data_path <- Tcas_path
  out_path_vst <- Tcas_path_out
  out_path_lfc <- Tcas_path_out_lfc
  metadata <- Tcas_metadata
}

#### read data for DE
## count data
count_data <- read.csv(data_path, header = TRUE, sep = "\t",comment.char = "#")
count_data <- count_data %>% select(-one_of('Chr','Start','End','Strand','Length')) 
# make the Geneid the row id 
count_data <- count_data %>% remove_rownames %>% column_to_rownames(var="Geneid")
## meta data
DE_metadata_all <- read.csv(metadata, header = TRUE, sep = ",",comment.char = "#")
# only include metadata when the samples are in the counts
meta_data <- DE_metadata_all[which((DE_metadata_all$Run %in% names(count_data))==TRUE), ]
meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var="Run")
# reorder so it's in the same order as the count data columns
meta_data <- meta_data[names(count_data),]
# check assumptions
all(colnames(count_data) %in% rownames(meta_data))


#### set up DESeq data structures
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = meta_data,
                              design = organ ~ sex)
# filter any counts less than 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# get normalized counts
dds <- estimateSizeFactors(dds)
# for the output files we might want this, but the new data structure doesn't work with DESeq2 any more
if (TRUE){
  # dds_vst <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
  # normalized_counts <- counts(dds_vst, normalized = TRUE)
  normalized_counts <- counts(dds, normalized = TRUE)
  write.table(normalized_counts, file=out_path_vst, sep = '\t', quote=F, col.names = NA)
}


#### run differential expression
dds <- DESeq(dds)
res <- results(dds)
summary(res)
write.table(res, file=out_path_lfc, sep = '\t', quote=F, col.names = NA)
## plot results
rld <- rlog(dds)
plotPCA(rld, intgroup="sex")







