#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -n 5
#SBATCH -p core
#SBATCH -t 1-00:00:00
#SBATCH -J blastp_for_synteny
#SBATCH -o blastp_for_synteny.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools blast/2.15.0+

## remove the species name prefixes I added for orthofinder
# sed 's/>B_siliquastri__B_siliquastri__B_siliquastri_/>/g' /proj/naiss2023-6-65/Milena/chapter2/protein_data/B_siliquastri.faa > /proj/naiss2023-6-65/Milena/chapter2/protein_data/B_siliquastri_original_header.faa
# sed 's/A_obtectus__A_obtectus__/>/g' /proj/naiss2023-6-65/Milena/chapter2/protein_data/A_obtectus.faa > /proj/naiss2023-6-65/Milena/chapter2/protein_data/A_obtectus_original_header.faa
# sed 's/>C_chinensis__C_chinensis__C_chinensis_/>/g' /proj/naiss2023-6-65/Milena/chapter2/protein_data/C_chinensis.faa > /proj/naiss2023-6-65/Milena/chapter2/protein_data/C_chinensis_original_header.faa
# sed 's/>Cmac_Lome_diverse_/>/g' /proj/naiss2023-6-65/Milena/chapter2/protein_data/C_maculatus_superscaffolded.faa > /proj/naiss2023-6-65/Milena/chapter2/protein_data/C_maculatus_superscaffolded_original_header.faa
# sed 's/>T_castaneum__T_castaneum__T_castaneum_/>/g' /proj/naiss2023-6-65/Milena/chapter2/protein_data/T_castaneum.faa > /proj/naiss2023-6-65/Milena/chapter2/protein_data/T_castaneum_original_header.faa

A_obtectus_proteins=/proj/naiss2023-6-65/Milena/chapter3/protein_data/A_obtectus_original_header.faa
B_siliquastri_proteins=/proj/naiss2023-6-65/Milena/chapter3/protein_data/B_siliquastri_original_header.faa
C_chinensis_proteins=/proj/naiss2023-6-65/Milena/chapter3/protein_data/C_chinensis_original_header.faa
C_maculatus_proteins=/proj/naiss2023-6-65/Milena/chapter3/protein_data/C_maculatus_superscaffolded_original_header.faa
T_castaneum_proteins=/proj/naiss2023-6-65/Milena/chapter3/protein_data/T_castaneum_original_header.faa
T_freemani_proteins=/proj/naiss2023-6-65/Milena/chapter3/protein_data/T_freemani_original_header.faa
C_magnifica_proteins=/proj/naiss2023-6-65/Milena/chapter3/protein_data/C_magnifica_original_header.faa

## --> re-run for new proteinfiles!

## make databases
for SPECIES1 in $T_freemani_proteins $C_magnifica_proteins # $A_obtectus_proteins $B_siliquastri_proteins $C_chinensis_proteins $C_maculatus_proteins $T_castaneum_proteins
do
    makeblastdb -in $SPECIES1 -dbtype prot
    echo " ---> done database ${SPECIES1}"
done

## -->

for SPECIES1 in $A_obtectus_proteins $B_siliquastri_proteins $C_chinensis_proteins $C_maculatus_proteins $T_castaneum_proteins $T_freemani_proteins $C_magnifica_proteins
do  

    SPECIES1_name="${SPECIES1##*/}"
    SPECIES1_name="${SPECIES1_name%.*}"

    for SPECIES2 in $A_obtectus_proteins $B_siliquastri_proteins $C_chinensis_proteins $C_maculatus_proteins $T_castaneum_proteins $T_freemani_proteins $C_magnifica_proteins
    do

        SPECIES2_name="${SPECIES2##*/}"
        SPECIES2_name="${SPECIES2_name%.*}"

        # if [[ "${SPECIES1_name}" == "${SPECIES2_name}" ]]
        # then
        #     continue
        # fi

        OUT_1v2="${SPECIES1_name}_vs_${SPECIES2_name}.blast"
        OUT_2v1="${SPECIES2_name}_vs_${SPECIES1_name}.blast"

        # the documentation says outfmt6 but I think they mean 8
        echo "blastp -query $SPECIES1 -db $SPECIES2 -out $OUT_1v2 -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6"
        echo " =========> ${OUT_1v2} done!"

        # # reverse already happens automatically in the nested for loop no need to implement explicitly
        # blastp -query $SPECIES2 -db $SPECIES1 -out $OUT_2v1 -num_threads 5 -num_alignments 5 -evalue 1e-10  -outfmt 6
        # echo " =========> ${OUT_2v1} done!"

    done
done