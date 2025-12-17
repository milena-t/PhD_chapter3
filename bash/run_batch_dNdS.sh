#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -J run_dNdS
#SBATCH -o run_dNdS.log


IN_LIST=$@ # space separated list of filenames, do like 100 at a time

CHR_TYPE=X

RUN_DIR=/proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/src/blast_BRH/
if [ -d "$RUN_DIR" ]; then
    module load FastTree/2.2-GCCcore-13.3.0
    IN_DIR=/proj/naiss2023-6-65/Milena/chapter3/dNdS_calculations/brh_sequences_${CHR_TYPE}/
    PAL2NAL=/proj/naiss2023-6-65/Lila/beetle_genomes/pal2nal.v14/pal2nal.pl
    CLUSTALO=/proj/naiss2023-6-65/Milena/software_install/clustal_omega/clustal-omega-1.2.4/bin/clustalo
    PAML=/sw/bioinfo/paml/4.10.7/rackham/bin/codeml
    FASTTREE=FastTree
else
    RUN_DIR=/Users/miltr339/work/PhD_code/PhD_chapter3/src/blast_BRH/
    IN_DIR=/Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_sequences_${CHR_TYPE}/
    PAL2NAL=/Users/miltr339/work/pal2nal.v14/pal2nal.pl
    CLUSTALO=/Users/miltr339/work/clustal-omega-1.2.4/src/clustalo
    PAML=/Users/miltr339/work/paml/src/codeml
    FASTTREE=/Users/miltr339/Desktop/FastTree
fi

for FILE in $IN_LIST # D_carinulata_D_sublineata_X-linked_ortholog_2.fasta #  
do
    # bash run_dNdS_pelle.sh $FILE

    python3 "${RUN_DIR}calculate_pairwise_dNdS.py" \
        --cds "${IN_DIR}${FILE}" \
        --pal2nalbin $PAL2NAL \
        --codeml \
        --codemlbin $PAML \
        --clustalbin $CLUSTALO \
        --fasttreebin $FASTTREE \
        --codemlbin $PAML \
        --overwrite \
        --verbose
        
    echo ">>> DONE ${FILE}"

done