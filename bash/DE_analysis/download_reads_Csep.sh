#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -n 8
#SBATCH -t 1:00:00
#SBATCH -J download_sequences
#SBATCH -o download_sequences.log
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load sratools/3.0.7
module load pigz

# Define directories
BASE=/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella
SRA=$BASE/raw_data/sra
FASTQ=$BASE/raw_data/fastq

mkdir -p $SRA $FASTQ

prefetch $1 --output-directory $SRA # Download .sra file
fasterq-dump $SRA/$1/$1.sra -O $FASTQ --split-files --threads 8 # Convert to FASTQ (paired-end, gzipped)

gunzip "${FASTQ}"/*