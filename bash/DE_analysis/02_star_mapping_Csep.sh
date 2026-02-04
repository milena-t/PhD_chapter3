#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH -t 20:00:00
#SBATCH -J STAR_mapping_Csep
#SBATCH -o STAR_mapping_Csep.log
#SBATCH --mail-type=ALL

module load  STAR/2.7.11b-GCC-13.3.0

GENOME=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_septempunctata/assembly_genomic.fna.masked
INDEX=/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella/star_index
ANNOT=/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella/GCF_907165205.1_icCocSept1.1_genomic.gtf
OUT_DIR=/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella/star_mapping
RNA_DIR=/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella/raw_data/trimmed_fastp

TCAS_GENOME=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/T_castaneum/assembly_genomic.fna.masked


#Create the STAR genome index
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir "$INDEX" \
     --genomeFastaFiles "$GENOME" \
     --sjdbGTFfile "$ANNOT" \
     --sjdbOverhang 149 #max read length 150-1

#Align reads and produce gene counts
for R1 in "$RNA_DIR"/*_1_trimmed.fastq.gz; do
  R2="${R1/_1_trimmed.fastq.gz/_2_trimmed.fastq.gz}"
  echo $R1
  echo $R2
  SAMPLE="$(basename "$R1" "_1_trimmed.fastq.gz")"
  STAR --runThreadN 16 \
       --genomeDir "$INDEX" \
       --readFilesIn "$R1" "$R2" \
       --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate \
       --quantMode GeneCounts \
       --twopassMode Basic \
       --outFilterMultimapNmax 20 \
       --limitBAMsortRAM 20000000000 \
       --outFileNamePrefix "$OUT_DIR/${SAMPLE}_"
done