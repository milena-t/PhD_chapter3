#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 8
#SBATCH --mem=72G
#SBATCH -t 3-00:00:00
#SBATCH -J align_outgroup
#SBATCH -o align_outgroup.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se


LAST_PATH=/gorilla/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/software_install/last-aligner/last-1651/bin/
ASS_DIR=/gorilla/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/revision/MK_test/npstat_MK_test/assemblies
SPECIES=B_siliquastri
cd $ASS_DIR

# make reference database and train parameters on query species
# ${LAST_PATH}lastdb -P8 Cmacdb C_maculatus.fasta.masked
${LAST_PATH}last-train -P8 --revsym -C2 Cmacdb ${SPECIES}.fasta.masked > C_maculatus_${SPECIES}.train
# do alignment
${LAST_PATH}lastal -P8 -m2 -D1e9 -C2 --split -p C_maculatus_${SPECIES}.train Cmacdb ${SPECIES}.fasta.masked | ${LAST_PATH}last-split -r | ${LAST_PATH}maf-linked - > out_${SPECIES}.maf

${LAST_PATH}last-dotplot out_${SPECIES}.maf out_${SPECIES}.png