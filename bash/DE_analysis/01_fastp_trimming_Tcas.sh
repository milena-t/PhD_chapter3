#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH -t 2:00:00
#SBATCH -J trimming_fastp_Tcas
#SBATCH -o trimming_fastp_Tcas.log
#SBATCH --mail-type=ALL

module load fastp/1.0.1-GCC-13.3.0 MultiQC/1.28-foss-2024a

# Directory with reads
reads_dir="/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Tribolium/raw_data/fastq"

# Output directory for results
output_dir="/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Tribolium/raw_data/trimmed_fastp"
output_reports="/proj/naiss2023-6-65/Milena/chapter3/RNAseq/Tribolium/raw_data/trimmed_reports"


for SRR_num in  SRR1647922 SRR1647927 SRR1647933 SRR1647925 SRR1647926 SRR1647936
do
    # Infer the R1 and R2 filename
    r1="${reads_dir}/${SRR_num}_1.fastq"
    r2="${reads_dir}/${SRR_num}_2.fastq"

    # Defining output files
    out_r1="${output_dir}/${SRR_num}_1_trimmed.fastq.gz"
    out_r2="${output_dir}/${SRR_num}_2_trimmed.fastq.gz"

    echo "Running fastp on $SRR_num ..."
    fastp -w 16 -i "$r1" -I "$r2" -o "$out_r1" -O "$out_r2" -h "$output_reports/${SRR_num}fastp.html"
done

multiqc $output_reports