#!/bin/bash -l

BLAST_OUTDIR_UPPMAX=/proj/naiss2023-6-65/Milena/chapter3/all_vs_all_blastp
BLAST_OUTDIR=/Users/miltr339/work/pairwise_blast_chapter_2_3
ANN_DIR=/Users/miltr339/work/chapter3/isoform_filtered_native_annotations
SCRIPT_PATH=/Users/miltr339/work/PhD_code/PhD_chapter3/src/blast_BRH

python3 $SCRIPT_PATH/get_blast_BRH.py \
--blast1 $BLAST_OUTDIR/D_sublineata_original_header_vs_D_carinulata_original_header.blast \
--blast2 $BLAST_OUTDIR/D_carinulata_original_header_vs_D_sublineata_original_header.blast \
--annotation1 $ANN_DIR/D_sublineata.gff \
--annotation2 $ANN_DIR/D_carinulata.gff \
--X_contigs1 NC_079485.1 \
--X_contigs2 NC_079460.1 \
--Y_contigs1 NC_079486.1 \
--Y_contigs2 NC_079473.1