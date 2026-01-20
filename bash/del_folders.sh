#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -J del_old_Fasta
#SBATCH -o del_old_Fasta.log


rm -r /proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_sequences_A
rm -r /proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_sequences_X
