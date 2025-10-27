#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00
#SBATCH -J bedfiles_for_MCScanX
#SBATCH -o bedfiles_for_MCScanX.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools biopython/1.80-py3.10.8

# python3 /proj/snic2021-6-30/Milena/chapter2/PhD_chapter2/src/make_bedfile_for_MCScanX.py /proj/naiss2023-6-65/Milena/chapter2/native_annotations/A_obtectus.gff ao
# echo "done A_obtectus.gff"
# python3 /proj/snic2021-6-30/Milena/chapter2/PhD_chapter2/src/make_bedfile_for_MCScanX.py /proj/naiss2023-6-65/Milena/chapter2/native_annotations/B_siliquastri.gff bs
# echo "done B_siliquastri.gff"
# python3 /proj/snic2021-6-30/Milena/chapter2/PhD_chapter2/src/make_bedfile_for_MCScanX.py /proj/naiss2023-6-65/Milena/chapter2/native_annotations/C_chinensis.gff cc
# echo "done C_chinensis.gff"
# python3 /proj/snic2021-6-30/Milena/chapter2/PhD_chapter2/src/make_bedfile_for_MCScanX.py /proj/naiss2023-6-65/Milena/chapter2/native_annotations/C_maculatus_superscaffolded.gff cm
# echo "done C_maculatus_superscaffolded.gff"
# python3 /proj/snic2021-6-30/Milena/chapter2/PhD_chapter2/src/make_bedfile_for_MCScanX.py /proj/naiss2023-6-65/Milena/chapter2/native_annotations/T_castaneum.gff tc
# echo "done T_castaneum.gff"
# python3 /proj/snic2021-6-30/Milena/chapter2/PhD_chapter2/src/make_bedfile_for_MCScanX.py /proj/naiss2023-6-65/Milena/chapter3/species_annotations/C_magnifica/braker_isoform_filtered.gff mc #unsure if it gets fucked if i have cm twice
# echo "done C_magnifica.gff"
# python3 /proj/snic2021-6-30/Milena/chapter2/PhD_chapter2/src/make_bedfile_for_MCScanX.py /proj/naiss2023-6-65/Milena/chapter3/species_annotations/T_freemani/braker_isoform_filtered.gff tf
# echo "done C_freemani.gff"
python3 /proj/snic2021-6-30/Milena/chapter2/PhD_chapter2/src/make_bedfile_for_MCScanX.py /proj/naiss2023-6-65/Milena/chapter2/native_annotations/C_septempunctata.gff cs
echo "done C_septempunctata.gff"