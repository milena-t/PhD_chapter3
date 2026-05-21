#!/bin/sh

TRANSEQ_PATH=/proj/naiss2023-6-65/Milena/software_install/emboss/EMBOSS-6.6.0/EMBOSS-6.6.0/bin/transeq
PAL2NAL_UPPMAX=/proj/naiss2023-6-65/Lila/beetle_genomes/pal2nal.v14/pal2nal.pl
CLUSTALO_UPPMAX=/proj/naiss2023-6-65/Milena/software_install/clustal_omega/clustal-omega-1.2.4/bin/clustalo
PAML_UPPMAX=/sw/bioinfo/paml/4.10.7/rackham/bin/codeml

CDS_dir=/Users/milena/work/pairwise_blast_chapter_2_3/brh_seqs_rev_test_A

cd $CDS_dir

for CDS_FASTA in ${CDS_dir}/*
do

python3 /Users/milena/work/PhD_code/PhD_chapter3/src/blast_BRH/calculate_pairwise_dNdS.py \
    --cds $CDS_FASTA \
    --pal2nalbin /Users/milena/work/pal2nal.v14/pal2nal.pl \
    --codeml \
    --branch_pairwise \
    --codemlbin /Users/milena/work/paml/src/codeml \
    --clustalbin /Users/milena/work/clustal-omega-1.2.4/src/clustalo \
    --fasttreebin /Users/milena/Desktop/FastTree \
    --verbose \
    --overwrite
    break
done
    # --yn00 \
    # --yn00bin /Users/miltr339/work/paml/src/yn00 \

