#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -J run_dNdS
#SBATCH -o run_dNdS.log


rm -r /proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_sequences_A
rm -r /proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_sequences_X
