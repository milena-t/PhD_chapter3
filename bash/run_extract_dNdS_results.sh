#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -J run_extract_dNdS_results
#SBATCH -o run_extract_dNdS_results.log

module load Biopython/1.84-foss-2024a

echo "python3 ../src/blast_BRH/extract_dNdS_results.py" 
echo 
python3 ../src/blast_BRH/extract_dNdS_results.py