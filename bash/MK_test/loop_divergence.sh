#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=52G
#SBATCH -t 3:00:00
#SBATCH -J sp_aln_loop
#SBATCH -o sp_aln_loop.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load SAMtools/1.22.1-GCC-13.3.0 bwa-mem2/2.3-GCC-13.3.0

# Directory with reads
reads_dir="/proj/coleoptera-genomics-2025/snic2021-6-30/Martyna/PoolSeq/reads_trimmed"
mapped_dir="/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/revision/MK_test/bwa_mapping"

cd $mapped_dir

### index reference before mapping 
# CMAC_index=C_maculatus_superscaffolded_index_for_bwa
# ref_dir=/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/annotation_pipeline/only_orthodb_annotation/C_maculatus_superscaffolded
# bwa-mem2 index -p $CMAC_index "${ref_dir}/assembly_genomic.fna.masked"

for SPECIES in A_obtectus B_siliquastri C_chinensis C_magnifica C_septempunctata T_castaneum T_freemani
do 
    sbatch -J "${SPECIES}_divergence_vcf" -o "${SPECIES}_divergence_vcf.out" -t 3-00:00:00 /proj/coleoptera-genomics-2025/snic2021-6-30/Milena/chapter3/PhD_chapter3/bash/MK_test/divergence_vcf.sh $SPECIES
done
 