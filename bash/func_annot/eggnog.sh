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
    # Limit annotation to Arthropods (6656)
    # python3 ${INDIR}eggnog/create_dbs.py -m diamond -y --taxids 6656 --dbname arthropoda --data_dir ${INDIR}eggnog/eggnog_db

# they say to use this one but it is is outdated, they use an old url in the download script that is not found
    # python3 ${INDIR}eggnog/download_eggnog_data.py -y -H -d 6656 --data_dir ${INDIR}eggnog/eggnog_db
# use these instead, they are the same as the script but the correct url, run in OutDIR
    # cd $OutDIR
    # wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz && echo Decompressing... && gunzip eggnog.db.gz
    # wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.taxa.tar.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz && echo Decompressing... && tar -zxf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz
    # wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog_proteins.dmnd.gz http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz && echo Decompressing... && gunzip eggnog_proteins.dmnd.gz
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
