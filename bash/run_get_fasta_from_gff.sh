#!/bin/bash -l

# cd /proj/naiss2023-6-65/Milena/coleoptera_sequences/d_sublineata/
# sbatch -J d_sublineata_get_protein_seqs.log /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/get_fasta_from_gff.sh GCF_026230105.1_icDioSubl1.1_genomic_isoform_filtered.gff GCF_026230105.1_icDioSubl1.1_genomic.fna

# cd /proj/naiss2023-6-65/Milena/coleoptera_sequences/d_carinulata/
# sbatch -J d_carinulata_get_protein_seqs.log /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/get_fasta_from_gff.sh GCF_026250575.1_icDioCari1.1_genomic_isoform_filtered.gff GCF_026250575.1_icDioCari1.1_genomic.fna

cd /proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_Lome_diverse/
sbatch -J Cmac_superscaffolded_get_protein_seqs.log /proj/naiss2023-6-65/Milena/chapter3/PhD_chapter3/bash/get_fasta_from_gff.sh /proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_Lome_diverse/Cmac_Lome_diverse/braker/braker_isoform_filtered.gff assembly_genomic.fna.masked