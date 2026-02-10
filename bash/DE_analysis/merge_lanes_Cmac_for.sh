#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J merge_forw_lanes_Cmac
#SBATCH -o merge_forw_lanes_Cmac.log
#SBATCH --mail-type=ALL

cd /proj/naiss2023-6-65/Milena/chapter3/RNAseq/maculatus/raw_data/trimmed_fastp

## forward reads
zcat SRR3113348_1_trimmed.fastq.gz SRR3113384_1_trimmed.fastq.gz | gzip > lanes_merged_AVf1_1_trimmed.fastq.gz
zcat SRR3113350_1_trimmed.fastq.gz SRR3113385_1_trimmed.fastq.gz | gzip > lanes_merged_AVf2_1_trimmed.fastq.gz
zcat SRR3113352_1_trimmed.fastq.gz SRR3113386_1_trimmed.fastq.gz | gzip > lanes_merged_AVf3_1_trimmed.fastq.gz
zcat SRR3113420_1_trimmed.fastq.gz SRR3113365_1_trimmed.fastq.gz | gzip > lanes_merged_AVm1_1_trimmed.fastq.gz
zcat SRR3113421_1_trimmed.fastq.gz SRR3113366_1_trimmed.fastq.gz | gzip > lanes_merged_AVm2_1_trimmed.fastq.gz
zcat SRR3113422_1_trimmed.fastq.gz SRR3113367_1_trimmed.fastq.gz | gzip > lanes_merged_AVm3_1_trimmed.fastq.gz
zcat SRR3113345_1_trimmed.fastq.gz SRR3113380_1_trimmed.fastq.gz | gzip > lanes_merged_HtVf1_1_trimmed.fastq.gz
zcat SRR3113346_1_trimmed.fastq.gz SRR3113382_1_trimmed.fastq.gz | gzip > lanes_merged_HtVf2_1_trimmed.fastq.gz
zcat SRR3113347_1_trimmed.fastq.gz SRR3113383_1_trimmed.fastq.gz | gzip > lanes_merged_HtVf3_1_trimmed.fastq.gz
zcat SRR3113364_1_trimmed.fastq.gz SRR3113361_1_trimmed.fastq.gz | gzip > lanes_merged_HtVm1_1_trimmed.fastq.gz
zcat SRR3113381_1_trimmed.fastq.gz SRR3113362_1_trimmed.fastq.gz | gzip > lanes_merged_HtVm2_1_trimmed.fastq.gz
zcat SRR3113418_1_trimmed.fastq.gz SRR3113363_1_trimmed.fastq.gz | gzip > lanes_merged_HtVm3_1_trimmed.fastq.gz
echo ">>>>>>>> forward reads done!" 
