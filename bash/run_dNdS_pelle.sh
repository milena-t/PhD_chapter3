#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -J run_dNdS
#SBATCH -o run_dNdS.log

# ml FastTree/2.2-GCCcore-13.3.0 Biopython/1.84-foss-2024a PAML/4.10.9-GCCcore-13.3.0
INFILE=$1

PAL2NAL_UPPMAX=/proj/naiss2023-6-65/Lila/beetle_genomes/pal2nal.v14/pal2nal.pl
CLUSTALO_UPPMAX=/proj/naiss2023-6-65/Milena/software_install/clustal_omega/clustal-omega-1.2.4/bin/clustalo
PAML_UPPMAX=/sw/bioinfo/paml/4.10.7/rackham/bin/codeml
PAML_CONFIG=/sw/bioinfo/paml/4.10.7/rackham/examples/codeml.ctl
IN_DIR=/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_sequences_X/

python3 /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/src/blast_BRH/calculate_pairwise_dNdS.py \
    --cds "${IN_DIR}${INFILE}" \
    --pal2nalbin $PAL2NAL_UPPMAX \
    --codeml \
    --branch_model \
    --config_path $PAML_CONFIG \
    --clustalbin $CLUSTALO_UPPMAX \
    
    # --codemlbin $PAML_UPPMAX \
    # --verbose \
    # --overwrite
    # --fasttreebin /Users/miltr339/Desktop/FastTree \
    # --yn00 \
    # --yn00bin /Users/miltr339/work/paml/src/yn00 \

