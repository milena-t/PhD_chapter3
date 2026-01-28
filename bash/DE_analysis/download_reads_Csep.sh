#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 8
#SBATCH -t 1:00:00
#SBATCH -J download_sequences
#SBATCH -o download_sequences_Csep.log
#SBATCH --mail-type=ALL

module load pigz/2.8-GCCcore-13.3.0

# Define directories
BASE=/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella
SRA=$BASE/raw_data/sra
FASTQ=$BASE/raw_data/fastq
SRAPATH=/proj/naiss2023-6-65/Milena/software_install/sra_tools/sratoolkit.3.3.0-ubuntu64/bin/

mkdir -p $SRA $FASTQ

${SRAPATH}prefetch $1 --output-directory $SRA # Download .sra file
${SRAPATH}fasterq-dump $SRA/$1/$1.sra -O $FASTQ --split-files --threads 8 # Convert to FASTQ (paired-end, gzipped)

gunzip "${FASTQ}"/*