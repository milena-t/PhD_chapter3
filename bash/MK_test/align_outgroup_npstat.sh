#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 8
#SBATCH --mem=72G
#SBATCH -t 30:00:00
#SBATCH -J align_outgroup
#SBATCH -o align_outgroup.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

ml GCCcore/14.3.0 BEDTools/2.31.1-GCC-13.3.0

LAST_PATH=/gorilla/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/software_install/last-aligner/last-1651/bin/
ASS_DIR=/gorilla/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/revision/MK_test/npstat_MK_test/assemblies
SPECIES=$1
SPECIES=C_chinensis
cd $ASS_DIR

if [ 1 == 1 ] ; then
    
    # make reference database and train parameters on query species
    # ${LAST_PATH}lastdb -P8 Cmacdb C_maculatus.fasta.masked
    
    echo "running ------- last-train "
    # ${LAST_PATH}last-train -P8 --revsym -C2 Cmacdb ${SPECIES}.fasta.masked > C_maculatus_${SPECIES}.train
    echo "done ---------- last-train "

    # do alignment
    echo "running ------- lastal "
    ${LAST_PATH}lastal -P8 -m2 -D1e9 -C2 --split -p C_maculatus_${SPECIES}.train Cmacdb ${SPECIES}.fasta.masked | ${LAST_PATH}last-split -r | ${LAST_PATH}maf-linked - > out_${SPECIES}.maf
    echo "done ---------- lastal "

    echo "running ------- last-dotplot "
    ${LAST_PATH}last-dotplot out_${SPECIES}.maf out_${SPECIES}.png
    echo "done ---------- last-dotplot "
fi

echo "running ------- check coverage"
${LAST_PATH}maf-convert bed out_${SPECIES}.maf > aln_${SPECIES}.bed
sort -k1,1 -k2,2n aln_${SPECIES}.bed > aln_${SPECIES}.sorted.bed
bedtools merge -i aln_${SPECIES}.sorted.bed > aln_${SPECIES}.merged.bed
awk '{sum += $3-$2} END {print sum}' aln_${SPECIES}.merged.bed