#!/bin/sh

TFRE_PATH=/Users/miltr339/work/repeatmasking/Tfre/

python3 /Users/miltr339/work/PhD_code/ReVis/src/ReVis/ReVis.py \
    --masker_outfile "${TFRE_PATH}Tfre_GCA_022388455.1.fasta.ori.out" \
    --masker_out_gff "${TFRE_PATH}Tfre_GCA_022388455.1.fasta.out.gff" \
    --out_dir $TFRE_PATH \
    --species_name T_freemani \
    --window_length 5e5 \
    --verbose \
    --plot
    #--plot_overlap_filtered \