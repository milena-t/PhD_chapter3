### Everything is mostly based off of this blog post:
# https://ashleyschwartz.com/posts/2023/05/deseq2-tutorial
# and i use vst for count normalization: https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/varianceStabilizingTransformation


# install DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2") # compiled from source 
# install tidyverse (ggplot2 is in tidyverse)
install.packages("tidyverse")
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

#### metadata relative paths
Cmac_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Cmac_full_SRR_list.csv'
Csep_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Csep_full_SRR_list.csv'
Tcas_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Tcas_full_SRR_list.csv'

#### read data for DE
if (TRUE){
  data_path <- Cmac_path
  out_path_vst <- Cmac_path_out
  metadata <- Cmac_metadata
}
if (FALSE){
  data_path <- Csep_path
  out_path_vst <- Csep_path_out
  metadata <- Csep_metadata
}
if (FALSE){
  data_path <- Tcas_path
  out_path_vst <- Tcas_path_out
  metadata <- Tcas_metadata
}
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
if (FALSE){
  dds_vst <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
  write.table(normalized_counts, file=out_path_vst, sep = '\t', quote=F, col.names = NA)
}


#### run differential expression
dds <- DESeq(dds)
res <- results(dds)
summary(res)






