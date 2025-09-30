#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00
#SBATCH -J convert_gbff_to_gff
#SBATCH -M rackham
#SBATCH -o convert_gbff_to_gff.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools biopython/1.80-py3.10.8

ANN_DIR=/proj/naiss2023-6-65/Milena/chapter3/species_annotations

python3 /proj/naiss2023-6-65/Milena/chapter3/convert_genbank_to_gff3.py \
-i ${ANN_DIR}/Cmag_GCA_965644565.1.gbff \
-o {ANN_DIR}/Cmag_GCA_965644565.1.gff \
--no_fasta

