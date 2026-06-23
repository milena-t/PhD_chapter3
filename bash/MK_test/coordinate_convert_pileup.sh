#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=72G
#SBATCH -t 3-00:00:00
#SBATCH -J liftover_bam
#SBATCH -o liftover_bam.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se


cd /gorilla/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/revision/MK_test/npstat_MK_test

CrossMap bam utg_to_superscaffold.chain UF-2981-ANC_S1_L001.bam UF-2981-ANC_S1_L001.superscaffolded.bam