#!/bin/bash -l

BLAST_OUTDIR_UPPMAX=/proj/naiss2023-6-65/Milena/chapter3/all_vs_all_blastp
BLAST_OUTDIR=/Users/miltr339/work/pairwise_blast_chapter_2_3
ANN_DIR=/Users/miltr339/work/chapter3/isoform_filtered_native_annotations

python3 src/get_blast_BRH.py \
--blast1 $BLAST_OUTDIR/A_obtectus_original_header_vs_B_siliquastri_original_header.blast \
--blast2 $BLAST_OUTDIR/B_siliquastri_original_header_vs_A_obtectus_original_header.blast \
--annotation1 $ANN_DIR/A_obtectus.gff \
--annotation2 $ANN_DIR/B_siliquastri.gff \
--X_contigs1 CAVLJG010000002.1,CAVLJG010003236.1,CAVLJG010003544.1,CAVLJG010000099.1,CAVLJG010000155.1,CAVLJG010000244.1,CAVLJG010000377.1,CAVLJG010000488.1 \
--X_contigs2 X \
--Y_contigs1 CAVLJG010000343.1,CAVLJG010002896.1,CAVLJG010000233.1,CAVLJG010000566.1,CAVLJG010000588.1 \
--Y_contigs2 Y 