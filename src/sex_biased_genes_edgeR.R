## edgeR
## user guide here: https://bioconductor.org/packages//release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# install all dependencies
if (FALSE){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2") # compiled from source 
  BiocManager::install("edgeR")
  # install tidyverse (ggplot2 is in tidyverse)
  install.packages("tidyverse") 
  install.packages("gtable") 
}
library(tidyverse)
library(ggplot2)
library(edgeR)
library(limma)
library(DESeq2)

######################
#### load data
######################
if (TRUE){
  setwd("/Users/miltr339/work/chapter3/DE_analysis")
  Cmac_path <- "Cmac_gene_counts_short_headers.txt"
  Csep_path <- "Csep_gene_counts_short_headers.txt"
  Tcas_path <- "Tcas_gene_counts_short_headers.txt"
  ## normalized counts outfiles
  Cmac_path_out <- "/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_gene_counts_edgeR_normalized.txt"
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
  data_path <- Cmac_path
  out_path_vst <- Cmac_path_out
  out_path_lfc <- Cmac_path_out_lfc
  metadata <- Cmac_metadata
}

######################
##### read data ffrom files and check if metadata matches
######################
if (TRUE){
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
}

#####################
#### read into EdgeR datastructures and do some preprocessing
#####################
if (TRUE){
  # make sample groups that reflect tissue (organ) and sex
  groups <- interaction(meta_data$organ, meta_data$sex, sep="_")
  # read data into data structure
  DE_counts <- DGEList(counts=count_data, group=groups)
  #dim(DE_counts)
  #[1] 14861    12
  
  ## filter out low count genes
  keep <- filterByExpr(DE_counts, group=groups)
  DE_filtered <- DE_counts[keep, , keep.lib.sizes=FALSE]
  DE_filtered$samples$lib.size <- colSums(DE_filtered$counts) # recompute library sizes
  # dim(DE_filtered)
  # [1] 11337    12
  
  ## normalize by library sizes
  DE_normalized <- calcNormFactors(DE_filtered)
}
## export normalized counts to file
if (TRUE){
  log_norm_counts <- cpm(DE_normalized,normalized.lib.sizes = TRUE,log = TRUE,prior.count = 1)
  write.table(log_norm_counts,file = out_path_vst,sep = "\t",quote = FALSE)# ,col.names = NA)
  
}

#####################
#### set up DE analysis
#####################
if (TRUE){
  MM <- model.matrix(~0 + groups) # suppress intercept with '0+'
  rownames(MM) <- colnames(DE_normalized)
  # colnames(MM) <- c("abdomen_female", "abdomen_male", "head_female","head_male") ### check if the order is right before assigning!!!
  colnames(MM) <- c("abdomen_female", "head_female", "abdomen_male","head_male") ### check if the order is right before assigning!!!
  
  # Estimates common, trended and gene-wise dispersion, allowing for a possible trend with average count size 
  DE_dispersion <- estimateDisp(DE_normalized, MM)
  plotBCV(DE_dispersion)
  
  # define contrasts
  my.contrasts <- makeContrasts(
    Sex_a =	abdomen_female - abdomen_male, 
    Sex_h =	head_female - head_male, 
    levels=MM)
  
  # fit linear model
  fit <- glmFit(DE_dispersion, MM)
}


#####################
#### analyze contrasts
#####################
if (TRUE){
  lrt_a <- glmLRT(fit, contrast=my.contrasts[,"Sex_a"])
  summary(de<-decideTestsDGE(lrt_a, p.value=0.05, lfc=1))#, adjust.method="fdr") 
  '1*abdomen_female -1*abdomen_male --> downregulated=male-biased
  Down                               3265
  NotSig                             5279
  Up                                 2793'
  lrt_h <- glmLRT(fit, contrast=my.contrasts[,"Sex_h"])
  summary(de<-decideTestsDGE(lrt_h, p.value=0.05, lfc=1))#, adjust.method="fdr") 
  '1*head_female -1*head_male
  Down                         1331
  NotSig                       9201
  Up                            805'
  
  ### try to reproduce above numbers from topTags()
  # abdomen (correct!)
  if (FALSE){
    SexbAbd <- topTags(lrt_a, n=Inf)
    SexbAbd_dat <-SexbAbd@.Data[[1]] 
    length(SexbAbd_dat$logFC) # 11337 -> full correct dataset
    SexbAbd_sig <- SexbAbd_dat[SexbAbd_dat$FDR < 0.05, ]
    length(SexbAbd_sig$logFC) # 9353
    
    
    SexbAbd <- topTags(lrt_a, n=Inf, p.value=0.05)	
    SexbAbd_dat <-SexbAbd@.Data[[1]] 
    length(SexbAbd_dat$logFC) # 9353
    SexbAbd_up <- SexbAbd_dat[SexbAbd_dat$logFC > 1, ] # 2793
    SexbAbd_do <- SexbAbd_dat[SexbAbd_dat$logFC < -1, ] # 3265
  }
  
  ## save LFC
  ## LFC is female-male !!
  SexbAbd <- topTags(lrt_a, n=Inf)	
  write.table(SexbAbd,file = "/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_sex_DE_abdomen_edgeR.txt",sep = "\t",quote = FALSE)# ,col.names = NA)
  SexbHT <- topTags(lrt_h, n=Inf)
  write.table(SexbHT,file = "/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/Cmac_sex_DE_head_thorax_edgeR.txt",sep = "\t",quote = FALSE)# ,col.names = NA)
}

if (TRUE){
  ## smear plots for head and abdomen
  ## LFC is female-male !!
  # Abdomen
  de1<-decideTestsDGE(lrt_a, p.value=0.001)
  detags1 <- rownames(DE_dispersion)[as.logical(de1)]	
  png("/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/edgeR_analysis/smear_abdomen_LFC.png", width = 1500, height = 1200, res = 220, bg = "transparent")
  plotSmear(lrt_a, de.tags=detags1, smearWidth=0.1, smooth.scatter=F, xlim=c(), main="abdomen LFC (female-male)")
  abline(h=c(-1, 1), col="blue")
  dev.off()
  # head+thorax
  de2<-decideTestsDGE(lrt_h, p.value=0.001)
  detags2 <- rownames(DE_dispersion)[as.logical(de2)]	
  png("/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/edgeR_analysis/smear_head_thorax_LFC.png", width = 1500, height = 1200, res = 220, bg = "transparent")
  plotSmear(lrt_h, de.tags=detags2, smearWidth=0.1, smooth.scatter=F, xlim=c(), main="head+thorax LFC (female-male)")
  abline(h=c(-1, 1), col="blue")
  dev.off()
}

### plot genes with significant sex bias
if(TRUE){
  sig_level = 0.001
  
  Res_a <- topTags(lrt_a, n=Inf)
  Res_h <- topTags(lrt_h, n=Inf)
  Res_aDE <- topTags(lrt_a, n=Inf, p.value=sig_level)
  Res_hDE <- topTags(lrt_h, n=Inf, p.value=sig_level)
  
  abdg <- Res_a$table	
  htg <- Res_h$table
  
  abdg$gene <- row.names(abdg)
  htg$gene <- row.names(htg)
  
  df <- merge(abdg, htg, by="gene", all=T)
  df <- subset(df, FDR.x<sig_level | FDR.y<sig_level )
  # dim(df)
  # [1] 8797   11
  
  # factor to show significance level in each tissue
  df$category <- ifelse((df$FDR.x <sig_level) & (df$FDR.y >sig_level), "AbdOnly", 
                        ifelse((df$FDR.x >sig_level) & (df$FDR.y < sig_level), "HTOnly", 
                               ifelse((df$FDR.x <sig_level) & (df$FDR.y <sig_level),"Both", 0)))

  abdLogFC <- df[,c("logFC.x")] # Take only the logFC
  htLogFC <- df[,c("logFC.y")]
  abdLogFC <- as.numeric(abdLogFC)
  htLogFC <- as.numeric(htLogFC)
  
  category <- df[,c("category")]
  newDF <- as.data.frame(cbind(abdLogFC, htLogFC))	
  newDF$category <- as.factor(category)
  
  ## enter a pact with the devil
  g = ggplot(data= newDF, aes(x= htLogFC, y= abdLogFC)) +
    geom_point(aes(colour=factor(category)), alpha=0.6, size=1.75) + 
    scale_color_manual(values = c("#F2933A", "#4d7298", "#b82946"),name =sprintf("Differential expression (q<%f)", sig_level), 
                       labels=c("Abdomen only", "Whole body", "Head+Thorax only")) +
    xlab("log2(female-male) head+thorax") + ylab("log2(female-male) abdomen") + 
    ggtitle("Sex-biased expression in reproductive and somatic tissues") + 
    geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed") + 
    geom_vline(aes(xintercept=0), colour="#990000", linetype="dashed") +  
    geom_hline(aes(yintercept=1), colour="darkgrey", linetype ="dashed") + 
    geom_hline(aes(yintercept=-1), colour="darkgrey", linetype ="dashed") + 
    geom_vline(aes(xintercept =1), colour="darkgrey", linetype ="dashed") + 
    geom_vline(aes(xintercept =-1), colour="darkgrey", linetype="dashed") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  png("/Users/miltr339/work/PhD_code/PhD_chapter3/data/DE_analysis/edgeR_analysis/DE_tissues.png", width = 1600, height = 1300, res = 250,bg = "transparent")
  g			
  dev.off()
  
  
}


