#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH -t 3:00:00
#SBATCH -J alignment
#SBATCH -o alignment.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

ml GCC/14.3.0

for SPECIES in A_obtectus B_siliquastri C_chinensis C_magnifica C_septempunctata T_castaneum T_freemani
do
    echo ""
    echo "------------ ${SPECIES}"
    sbatch -J "${SPECIES}_alignment" -o "${SPECIES}_alignment.out" -t 3-00:00:00 /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/PhD_chapter3/bash/MK_test/align_outgroup_npstat.sh $SAMPLE
done
 