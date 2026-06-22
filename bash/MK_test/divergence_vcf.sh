#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=72G
#SBATCH -t 30:00:00
#SBATCH -J divergence_vcf
#SBATCH -o divergence_vcf.log

ml minimap2/2.30-GCCcore-13.3.0 k8/1.2-GCCcore-13.3.0 tabixpp/1.1.2-GCC-13.3.0

ASS_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/revision/MK_test/divergence_vcf/assemblies/
PAFTOOLS_PATH=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/software_install/minimap2/

SPECIES=$1

echo "====================== ${SPECIES} ======================"
QUERY_ASSEMBLY="${ASS_DIR}${SPECIES}.fasta.masked"
REF_CMAC_ASSEMBLY="${ASS_DIR}C_maculatus.fasta.masked"

## mapping
# Cmac needs to be ref so that it uses the Cmac contig coordinates
minimap2 -c -x asm10 $REF_CMAC_ASSEMBLY $QUERY_ASSEMBLY > aln_${SPECIES}.paf

## sorting alignment
sort -k6,6 -k8,8n aln_${SPECIES}.paf > aln_${SPECIES}.srt.paf

## variant calling (needs k8)
# I got paftools with wget https://raw.githubusercontent.com/lh3/minimap2/refs/heads/master/misc/paftools.js which is 2.31 and not 2.30 like the module
${PAFTOOLS_PATH}paftools.js call -f ingroup_reference.fasta -L 5000 -l 5000 aln_${SPECIES}.srt.paf > ${SPECIES}_C_maculatus_divergence.vcf

## sorting vcf
bgzip ${SPECIES}_C_maculatus_divergence.vcf && tabix -p vcf ${SPECIES}_C_maculatus_divergence.vcf.gz
