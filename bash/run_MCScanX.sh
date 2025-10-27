#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 1
#SBATCH -p core
#SBATCH -t 1:00:00
#SBATCH -J run_MCScanX
#SBATCH -o run_MCScanX.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# Bioinfo tools not needed on pelle
module load bioinfo-tools 

## test
# /proj/naiss2023-6-65/Milena/chapter2/MCScanX-1_0_0/MCScanX /proj/naiss2023-6-65/Milena/chapter2/MCScanX_test
# ./MCScanX test

## runs really fast for five coleopteran species
./MCScanX all_species
