#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 2:00:00
#SBATCH -J BUSCO
#SBATCH -o busco.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools BUSCO/5.5.0

# use the arthropoda_odb10 lineage database

# make softlinks of all assemblies
# for ASSEMBLY in $(echo /proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/*/assembly_genomic.fna) ; do echo $ASSEMBLY ; SPECIES="${ASSEMBLY%/assembly_genomic.fna}" ; echo $SPECIES ; SPECIES="${SPECIES#/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/}" ; echo $SPECIES ; ln -s $ASSEMBLY "${SPECIES}_genomic_assembly.fasta" ; done

# the working directory is created from scratch, so don't do a full path like /proj/naiss...
# WD=busco/

cd /proj/naiss2023-6-65/Milena/chapter3/species_annotations/BUSCO

busco -i T_freemani_isoform_filtered_proteins.faa -m proteins -l $BUSCO_LINEAGE_SETS/arthropoda_odb10/ -o T_freemani_busco -c 20 -f 
busco -i C_magnifica_isoform_filtered_proteins.faa -m proteins -l $BUSCO_LINEAGE_SETS/arthropoda_odb10/ -o C_magnifica_busco -c 20 -f 
