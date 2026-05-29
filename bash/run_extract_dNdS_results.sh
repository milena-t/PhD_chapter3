#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -t 10:00:00
#SBATCH --mail-type=ALL
#SBATCH -J revisions_extract_results
#SBATCH -o revisions_extract_results.log

module load Biopython/1.84-gfbf-2024a

echo "python3 /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/PhD_chapter3/src/blast_BRH/extract_dNdS_Results.py" 
echo 
python3 /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/PhD_chapter3/src/blast_BRH/extract_dNdS_Results.py