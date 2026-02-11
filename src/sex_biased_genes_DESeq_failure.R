# install all dependencies
if (FALSE){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2") # compiled from source 
  # install tidyverse (ggplot2 is in tidyverse)
  install.packages("tidyverse") 
  install.packages("gtable") 
  install.packages("pcadapt") 
}
library("DESeq2")
library(tidyverse)
library(ggplot2)
library(pcadapt)

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
Cmac_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Cmac_SRR_metadata.csv'
Csep_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Csep_full_SRR_list.csv'
Tcas_metadata <- '/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/metadata/Tcas_full_SRR_list.csv'

#### do one species of the three
if (TRUE){
  data_path <- Cmac_path
  out_path_vst <- Cmac_path_out
  out_path_lfc <- Cmac_path_out_lfc
  metadata <- Cmac_metadata
}
if (FALSE){
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
meta_data <- DE_metadata_all[which((DE_metadata_all$ID %in% names(count_data))==TRUE), ]
meta_data <- meta_data %>% remove_rownames %>% column_to_rownames(var="ID")
# reorder so it's in the same order as the count data columns
meta_data <- meta_data[names(count_data),]
# check assumptions
all(colnames(count_data) %in% rownames(meta_data))


#### set up DESeq data structures
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = meta_data,
                              design = ~0+sex) # the 0 suppresses the intercept
# filter any counts less than 10
keep <- rowSums(counts(dds) >=2) >= 3 # at least two counts in at least three samples
dds <- dds[keep,]
# get normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
write.table(normalized_counts, file=out_path_vst, sep = '\t', quote=F, col.names = NA)

if (FALSE){
  dds_vst <- rlog(dds, fitType = "local") #variance stabilizing transormation
  
  ## test PCA plotting with different tests
  # built-in function for DEseq2:
  plotPCA(dds_vst, intgroup="sex")
  plotPCA(dds_vst, intgroup="organ") # by organ i mean tissue
   
  ## try manually 
  # this gives completely different results, from both the above plotPCA and also the python PCA based on the normalized_counts
  # ?????????
  pca <- prcomp(normalized_counts)#, center = TRUE)# , scale. = TRUE)
  # Extract first two PCs for the loadings
  x <- pca$rotation[, "PC1"]
  y <- pca$rotation[, "PC2"]
  # Percent variance explained
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
  percent_var <- round(100 * var_explained, 1)
  # Plot points
  plot(x, y,pch = 19,xlab = paste0("PC1 (", percent_var[1], "%)"),ylab = paste0("PC2 (", percent_var[2], "%)"))
  text(x, y,labels = colnames(normalized_counts),pos = 3)
  
  ## PCadapt
  # again a new 
  expr <- t(normalized_counts)
  # expr <- normalized_counts
  pcadapt_input <- read.pcadapt(expr, type = "pcadapt")
  pc <- pcadapt(input = pcadapt_input, K = 10)
  percent_var <- round(100 * pc$singular.values^2 / sum(pc$singular.values^2), 1)
  plot(pc$scores[,1],pc$scores[,2],pch = 19,xlab = paste0("PC1 (", percent_var[1], "%)"),ylab = paste0("PC2 (", percent_var[2], "%)"))
  text(pc$scores[,1],pc$scores[,2],labels = colnames(normalized_counts),pos = 3,cex = 0.7)
  }



#### run differential expression
dds_DE <- DESeq(dds)
res <- results(dds_DE, alpha=0.05)
summary(res)

write.table(res, file=out_path_lfc, sep = '\t', quote=F, col.names = NA)

### plot results
## PCA
rld_DE<- rlog(dds_DE) # also very similar to VST, see function documentation
plotPCA(rld_DE, intgroup="sex")
## MA plot
plotMA(res, ylim=c(-5,5))

