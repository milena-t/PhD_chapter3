#!/bin/sh

ORTHODB_ARTHROPODA=/proj/naiss2023-6-65/Milena/annotation_pipeline/annotation_protein_data/OrthoDB_Arthropoda_v11.fa

# T. freemani annotation
TFRE_DIR=/proj/naiss2023-6-65/Milena/chapter3/species_assemblies
cd $TFRE_DIR

sbatch --job-name="T_freemani_annotation" --output="T_freemani_annotation.out" -t 5-00:00:00 \
/proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh \
T_freemani ${TFRE_DIR}/Tfre_GCA_022388455.1.fasta.masked \
$ORTHODB_ARTHROPODA SRR15965976 SRR15965980 SRR15965983 SRR15965985SRR14070854 SRR14070855SRR14070871 SRR14070870 

# SRR15965976 SRR15965980 SRR15965983 SRR15965985 \ # larvae RNA different stages and population densities
# SRR14070854 SRR14070855 \ # male adult
# SRR14070871 SRR14070870 \ # female adult