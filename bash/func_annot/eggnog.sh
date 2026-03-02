#!/bin/bash -l

#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH -t 50:00:00
#SBATCH -J eggnog_funcannot
#SBATCH -o eggnog_funcannot.log


## Ingo says:
    # Run eggnog-mapper with diamond, running with iterations similar to mmseqs2 with the final iteration on ultra-sensitive settings which should be slightly better than
    # mmseqs2 at s=7.5 (Buchfink, et al. (2025). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature methods)

# Load modules
module load eggnog-mapper/2.1.12-foss-2024a

INDIR=/proj/naiss2023-6-65/Milena/chapter3/Cmac_func_annot/
InFAA=${INDIR}Cmac_superscaffolded_Lome_braker_proteins.faa
OutDIR=${INDIR}eggnog
GFF=${INDIR}Cmac_superscaffolded_Lome_braker.gtf


## download eggnog database
# I did:
#   wget http://eggnog6.embl.de/download/novel_fams-1.0.1/novel_fams.dmnd.gz
# they also suggest:
python3 ${INDIR}eggnog/download_eggnog_data.py --data_dir ${INDIR}eggnog/eggnog_db

#Use/Make scratch dir
export scratchDIR=${SNIC_TMP}/Cmac_lome_eggnog_diamond

if [ -d ${scratchDIR} ]; then
    echo "Scratch directory already exists: ${scratchDIR}"
else
    mkdir $scratchDIR
    echo "created directory in SNIC_TMP: $scratchDIR"
fi

#Limit annotation to Arthropods (6656)
emapper.py \
    -i $InFAA --output_dir $OutDIR -o C_mac_eggnog_diamond \
    --cpu 16 \
    --scratch_dir $scratchDIR \
    --data_dir ${INDIR}eggnog/eggnog_db \
    -m diamond --tax_scope 6656 \
    --decorate_gff ${GFF} \
	--sensmode ultra-sensitive --dmnd_iterate yes --matrix BLOSUM62
