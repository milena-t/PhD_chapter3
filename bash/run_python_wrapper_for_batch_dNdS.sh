#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -J run_wrapper_dNdS
#SBATCH -o run_wrapper_dNdS.log

# module load Biopython/1.84-foss-2024a
ml FastTree/2.2-GCCcore-13.3.0 Biopython/1.84-gfbf-2024a argtable/2.13-GCCcore-13.3.0 PAML/4.10.9-GCCcore-13.3.0

echo "python3 /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/PhD_chapter3/src/blast_BRH/calculate_batch_pw_dNdS.py"
echo 
python3 /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/PhD_chapter3/src/blast_BRH/calculate_batch_pw_dNdS.py "X"