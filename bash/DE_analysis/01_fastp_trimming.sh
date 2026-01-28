#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -n 16
#SBATCH -t 2:00:00
#SBATCH -J trimming_fastp
#SBATCH -o trimming_fastp.log
#SBATCH --mail-type=ALL

module load bioinfo-tools fastp/0.23.4 MultiQC/1.22.2

# Directory with reads
reads_dir="/proj/naiss2023-6-65/Milena/chapter2/Tribolium_poolseq/raw_data/fastq"

# Output directory for results
output_dir="/proj/naiss2023-6-65/Milena/chapter2/Tribolium_poolseq/raw_data/trimmed_fastp"
output_reports="/proj/naiss2023-6-65/Milena/chapter2/Tribolium_poolseq/raw_data/trimmed_reports"


for SRR_num in  SRR23732240 SRR23732241 SRR23732242 SRR23732243 SRR23732244 SRR23732245
do
    # Infer the R1 and R2 filename
    r1="${reads_dir}/${SRR_num}_1.fastq"
    r2="${reads_dir}/${SRR_num}_2.fastq"

    # Defining output files
    out_r1="${output_dir}/${SRR_num}_1_trimmed.fastq.gz"
    out_r2="${output_dir}/${SRR_num}_2_trimmed.fastq.gz"

    echo "Running fastp on $sample ..."
    # fastp -w 16 -i "$r1" -I "$r2" -o "$out_r1" -O "$out_r2" -h "$output_reports/${SRR_num}fastp.html"
done

multiqc $output_reports