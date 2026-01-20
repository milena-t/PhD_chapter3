#!/bin/bash -l

BLAST_OUTDIR_UPPMAX=/proj/naiss2023-6-65/Milena/chapter3/all_vs_all_blastp
BLAST_OUTDIR=/Users/miltr339/work/pairwise_blast_chapter_2_3
ANN_DIR=/Users/miltr339/work/chapter3/isoform_filtered_native_annotations
SCRIPT_PATH=/Users/miltr339/work/PhD_code/PhD_chapter3/src/blast_BRH


python3 $SCRIPT_PATH/get_blast_BRH.py \
--blast1 $BLAST_OUTDIR/T_castaneum_original_header_vs_T_freemani_original_header.blast \
--blast2 $BLAST_OUTDIR/T_freemani_original_header_vs_T_castaneum_original_header.blast \
--annotation1 $ANN_DIR/T_castaneum.gff \
--annotation2 $ANN_DIR/T_freemani.gff \
--X_contigs1 NC_087403.1 \
--X_contigs2 CM039461.1 \
# --Y_contigs1 none \
# --Y_contigs2 none


python3 $SCRIPT_PATH/get_blast_BRH.py \
--blast1 $BLAST_OUTDIR/C_magnifica_original_header_vs_C_septempunctata_original_header.blast \
--blast2 $BLAST_OUTDIR/C_septempunctata_original_header_vs_C_magnifica_original_header.blast \
--annotation1 $ANN_DIR/C_magnifica.gff \
--annotation2 $ANN_DIR/C_septempunctata.gff \
--X_contigs1 OZ286750.1 \
--X_contigs2 NC_058198.1 \
# --Y_contigs1 none \
# --Y_contigs2 none