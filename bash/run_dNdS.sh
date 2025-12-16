#!/bin/sh

TRANSEQ_PATH=/proj/naiss2023-6-65/Milena/software_install/emboss/EMBOSS-6.6.0/EMBOSS-6.6.0/bin/transeq
PAL2NAL_UPPMAX=/proj/naiss2023-6-65/Lila/beetle_genomes/pal2nal.v14/pal2nal.pl
CLUSTALO_UPPMAX=/proj/naiss2023-6-65/Milena/software_install/clustal_omega/clustal-omega-1.2.4/bin/clustalo
PAML_UPPMAX=/sw/bioinfo/paml/4.10.7/rackham/bin/codeml

python3 /Users/miltr339/work/PhD_code/PhD_chapter3/src/blast_BRH/calculate_pairwise_dNdS.py \
    --cds /Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_sequences_X/B_siliquastri_C_maculatus_X-linked_ortholog_166.fasta \
    --pal2nalbin /Users/miltr339/work/pal2nal.v14/pal2nal.pl \
    --yn00 \
    --yn00bin /Users/miltr339/work/paml/src/yn00 \
    --clustalbin /Users/miltr339/work/clustal-omega-1.2.4/src/clustalo \
    --fasttreebin /Users/miltr339/Desktop/FastTree \
    --verbose \
    --overwrite
    
# --codemlbin /Users/miltr339/work/paml/src/codeml \

