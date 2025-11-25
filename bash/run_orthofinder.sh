#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 2:00:00
#SBATCH -J orthofinder_for_fastX
#SBATCH -o orthofinder_for_fastX.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools OrthoFinder/2.5.2

# the output directory will be created by orthofinder and should NOT already exist!
OUT_DIR=/proj/naiss2023-6-65/Milena/chapter3/orthofinder

# orthofinder runs with protein sequences here! nucleotide sequences require the -d flag
IN_FASTA_DIR=/proj/naiss2023-6-65/Milena/chapter3/protein_data/orthofinder_dir

# remove it if it does exist, this overwrites all preexisting results
rm -r $OUT_DIR

orthofinder.py -t 20 -a 20 -f $IN_FASTA_DIR -o $OUT_DIR
