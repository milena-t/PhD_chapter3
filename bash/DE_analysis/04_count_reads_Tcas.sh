#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH -t 1-00:00:00
#SBATCH -J Tcas_count_reads
#SBATCH -o Tcas_count_reads
#SBATCH --mail-type=ALL

module load Subread/2.1.1-GCC-13.3.0

#Directories
ANNOTATION="/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Tribolium/GCF_031307605.1_icTriCast1.1_genomic.gtf"
INPUT_DIR="/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Tribolium/star_mapping/picard_marked_indexed"
OUT_DIR="/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Tribolium/gene_counts"

mkdir -p "$OUT_DIR"

# Mode 1: Standard counting (unique mappers only)
# No multimapping, at the gene level, used for differential expression analysis
echo "=== Mode 1: Standard (unique mappers only) ==="
featureCounts -T 16 \
  -a "$ANNOTATION" \
  -o "$OUT_DIR/gene_counts_standard.txt" \
  -p -B -C \
  -g "gene_id" \
  -t "exon" \
  -s 2 \
  "$INPUT_DIR"/*_marked_duplicates.bam
