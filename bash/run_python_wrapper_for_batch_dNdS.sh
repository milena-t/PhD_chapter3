#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -J run_wrapper_dNdS
#SBATCH -o run_wrapper_dNdS.log

module load Biopython/1.84-foss-2024a

echo "python3 /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/src/blast_BRH/calculate_batch_pw_dNdS.py"
echo 
python3 /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/src/blast_BRH/calculate_batch_pw_dNdS.py