#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH -t 20:00:00
#SBATCH -J STAR_mapping_Csep
#SBATCH -o STAR_mapping_Csep.log
#SBATCH --mail-type=ALL

module load  STAR/2.7.11b-GCC-13.3.0

CSEP_GENOME=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_septempunctata/assembly_genomic.fna.masked
CSEP_INDEX=/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella/star_index
CSEP_ANNOT=/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Coccinella/GCF_907165205.1_icCocSept1.1_genomic.gtf

TCAS_GENOME=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/T_castaneum/assembly_genomic.fna.masked


#Create the STAR genome index
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir "$CSEP_INDEX" \
     --genomeFastaFiles "$CSEP_GENOME" \
     --sjdbGTFfile "$CSEP_ANNOT" \
     --sjdbOverhang 149 #max read length 150-1

#Align reads and produce gene counts
for R1 in "$RNA_DIR"/*_1.fastq.gz; do
  R2="${R1/_1.fastq.gz/_2.fastq.gz}"
  SAMPLE="$(basename "$R1" "_1.fastq.gz")"
  STAR --runThreadN 16 \
       --genomeDir "$GENOME_DIR" \
       --readFilesIn "$R1" "$R2" \
       --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate \
       --quantMode GeneCounts \
       --twopassMode Basic \
       --outFilterMultimapNmax 20 \
       --limitBAMsortRAM 20000000000 \
       --outFileNamePrefix "$OUT_DIR/${SAMPLE}_"
done