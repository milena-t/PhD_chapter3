#!/bin/sh
#SBATCH -A naiss2024-5-135
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 20:00
#SBATCH -J single_exon_stats
#SBATCH -o single_exon_stats.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools python3

# all orthoDB annotations

for SPECIES in  A_obtectus B_siliquastri C_chinensis C_maculatus_superscaffolded C_magnifica C_septempunctata T_castaneum T_freemani
do
    python3 /Users/miltr339/work/PhD_code/PhD_chapter3/src/plotting/calculate_single_exon_stats.py \
    "/Users/miltr339/work/chapter3/isoform_filtered_native_annotations/${SPECIES}.gff" \
    True > "/Users/miltr339/work/chapter3/single_exon_stats/${SPECIES}_single_exon_stats_with_transcript_list.txt"
    echo "-----> done ${SPECIES}"
done