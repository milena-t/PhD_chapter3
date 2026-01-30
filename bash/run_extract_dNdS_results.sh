#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH -J dNdS_extract_results
#SBATCH -o dNdS_extract_results.log

module load Biopython/1.84-foss-2024a

echo "python3 /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/src/blast_BRH/extract_dNdS_Results.py" 
echo 
python3 /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/src/blast_BRH/extract_dNdS_Results.py