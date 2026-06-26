#!/bin/bash -l
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH -J run_branch_model_bruchini
#SBATCH -o run_branch_model_bruchini.log


IN_LIST=$@ # space separated list of filenames, do like 100 at a time

CHR_TYPE=A

RUN_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/PhD_chapter3/src/blast_BRH/
if [ -d "$RUN_DIR" ]; then
    module load FastTree/2.2-GCCcore-13.3.0
    module load Biopython/1.84-gfbf-2024a
    module load argtable/2.13-GCCcore-13.3.0
    IN_DIR=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/dNdS_calculations/brh_sequences_${CHR_TYPE}/
    PAL2NAL=/proj/coleoptera-genomics-2025/snic2021-6-30/Lila/beetle_genomes/pal2nal.v14/pal2nal.pl
    CLUSTALO=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/software_install/clustal_omega/clustal-omega-1.2.4/bin/clustalo
    CODEML=/sw/bioinfo/paml/4.10.7/rackham/bin/codeml
    FASTTREE=FastTree
    PAML_CONFIG=/sw/bioinfo/paml/4.10.7/rackham/examples/codeml.ctl
else
    RUN_DIR=/Users/miltr339/work/PhD_code/PhD_chapter3/src/blast_BRH/
    IN_DIR=/Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/brh_sequences_${CHR_TYPE}/
    PAL2NAL=/Users/miltr339/work/pal2nal.v14/pal2nal.pl
    CLUSTALO=/Users/miltr339/work/clustal-omega-1.2.4/src/clustalo
    CODEML=/Users/miltr339/work/paml/src/codeml
    FASTTREE=/Users/miltr339/Desktop/FastTree
    PAML_CONFIG=/Users/miltr339/work/paml/examples/codeml.ctl 
fi

cwd=$(pwd)
for FILE in $IN_LIST
# for FILE in /Users/miltr339/work/pairwise_blast_chapter_2_3/brh_tables/bruchini_seqs_revision_X/all_bruchini/*.fasta
do
    
    cd  $cwd

    python3 "${RUN_DIR}calculate_pairwise_dNdS.py" \
        --cds "${FILE}" \
        --pal2nalbin $PAL2NAL \
        --codeml \
        --branch_model \
        --codeml_config_path $PAML_CONFIG \
        --clustalbin $CLUSTALO \
        --codemlbin $CODEML \
        --fasttreebin $FASTTREE \
        --overwrite \
        --verbose

    echo ">>> DONE ${FILE}"

done
        # --verbose