#!/bin/bash -l

#SBATCH -A uppmax2026-1-8
#SBATCH -c 8
#SBATCH --mem=6G
#SBATCH -t 1-00:00:00
#SBATCH -J interproscan
#SBATCH -o interproscan.log

###!!! not available on pelle yet
ml InterProScan/5.62-94.0

InFAA=/proj/snic2019-35-58/water_strider/ingo/data/intermediate/nobackup/02.annotation/03.functional/01.interpro/00.input/braker_nostop.faa
OutDIR=/proj/snic2019-35-58/water_strider/ingo/data/intermediate/nobackup/02.annotation/03.functional/01.interpro/01.output/

interproscan.sh -i $InFAA -cpu 8 -t p -dp -pa --goterms --iprlookup -d $OutDIR
