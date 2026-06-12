#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH -t 1-00:00:00
#SBATCH -J mapping_bwa_sample
#SBATCH -o mapping_bwa_sample.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se


module load SAMtools/1.22.1-GCC-13.3.0 bwa-mem2/2.3-GCC-13.3.0

CMAC_index=C_maculatus_superscaffolded_index_for_bwa

# Directory with reads
reads_dir="/proj/coleoptera-genomics-2025/snic2021-6-30/Martyna/PoolSeq/reads_trimmed"
mapped_dir="/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/revision/MK_test/bwa_mapping"

cd $mapped_dir

SAMPLE=$1
r1="${reads_dir}/${SAMPLE}_R1.trim.fastq.gz"
r2="${reads_dir}/${SAMPLE}_R2.trim.fastq.gz"
echo "======================>> Running bwa-mem on $sample ..."
bwa-mem2 mem -t 20 -P $CMAC_index $r1 $r2 | samtools view -u | samtools sort -o "${mapped_dir}/${SAMPLE}.bam"
