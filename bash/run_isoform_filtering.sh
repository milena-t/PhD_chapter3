#!/bin/bash -l

cd /proj/naiss2023-6-65/Milena/coleoptera_sequences/d_carinulata/
sbatch /Users/miltr339/work/PhD_code/PhD_chapter3/bash/isoform_filter_gff.sh GCF_026250575.1_icDioCari1.1_genomic.gff

cd /proj/naiss2023-6-65/Milena/coleoptera_sequences/d_carinulata/
sbatch /Users/miltr339/work/PhD_code/PhD_chapter3/bash/isoform_filter_gff.sh GCF_026250575.1_icDioCari1.1_genomic.gff