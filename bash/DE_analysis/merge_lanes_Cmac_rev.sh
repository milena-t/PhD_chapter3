#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J merge_rev_lanes_Cmac
#SBATCH -o merge_rev_lanes_Cmac.log
#SBATCH --mail-type=ALL

cd /proj/naiss2023-6-65/Milena/chapter3/RNAseq/maculatus/raw_data/trimmed_fastp

## reverse reads
zcat SRR3113348_2_trimmed.fastq.gz SRR3113384_2_trimmed.fastq.gz | gzip > lanes_merged_AVf_2_trimmed.fastq.gz
zcat SRR3113350_2_trimmed.fastq.gz SRR3113385_2_trimmed.fastq.gz | gzip > lanes_merged_AVf_2_trimmed.fastq.gz
zcat SRR3113352_2_trimmed.fastq.gz SRR3113386_2_trimmed.fastq.gz | gzip > lanes_merged_AVf_2_trimmed.fastq.gz
zcat SRR3113420_2_trimmed.fastq.gz SRR3113365_2_trimmed.fastq.gz | gzip > lanes_merged_AVm_2_trimmed.fastq.gz
zcat SRR3113421_2_trimmed.fastq.gz SRR3113366_2_trimmed.fastq.gz | gzip > lanes_merged_AVm_2_trimmed.fastq.gz
zcat SRR3113422_2_trimmed.fastq.gz SRR3113367_2_trimmed.fastq.gz | gzip > lanes_merged_AVm_2_trimmed.fastq.gz
zcat SRR3113345_2_trimmed.fastq.gz SRR3113380_2_trimmed.fastq.gz | gzip > lanes_merged_HtVf_2_trimmed.fastq.gz
zcat SRR3113346_2_trimmed.fastq.gz SRR3113382_2_trimmed.fastq.gz | gzip > lanes_merged_HtVf_2_trimmed.fastq.gz
zcat SRR3113347_2_trimmed.fastq.gz SRR3113383_2_trimmed.fastq.gz | gzip > lanes_merged_HtVf_2_trimmed.fastq.gz
zcat SRR3113364_2_trimmed.fastq.gz SRR3113361_2_trimmed.fastq.gz | gzip > lanes_merged_HtVm_2_trimmed.fastq.gz
zcat SRR3113381_2_trimmed.fastq.gz SRR3113362_2_trimmed.fastq.gz | gzip > lanes_merged_HtVm_2_trimmed.fastq.gz
zcat SRR3113418_2_trimmed.fastq.gz SRR3113363_2_trimmed.fastq.gz | gzip > lanes_merged_HtVm_2_trimmed.fastq.gz
echo ">>>>>>>> reverse reads done!" 