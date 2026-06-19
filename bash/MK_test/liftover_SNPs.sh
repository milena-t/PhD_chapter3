#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 1
#SBATCH --mem=72G
#SBATCH -t 3:00:00
#SBATCH -J liftover_vcf_to_superscaffold
#SBATCH -o liftover_vcf_to_superscaffold.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

ml picard/3.4.0-Java-17

# java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary -R superscaffolded_Cmac.fasta

java -jar $EBROOTPICARD/picard.jar LiftoverVcf \
--INPUT freebayes_ANC_only.vcf \
--OUTPUT freebayes_superscaffolded.vcf \
--CHAIN utg_to_superscaffold.chain \
--REJECT superscaffolded_rejected_variants.vcf \
--REFERENCE_SEQUENCE superscaffolded_Cmac.fasta