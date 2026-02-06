#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH -t 2:00:00
#SBATCH -J trimming_fastp_Cmac
#SBATCH -o trimming_fastp_Cmac.log
#SBATCH --mail-type=ALL

module load fastp/1.0.1-GCC-13.3.0 MultiQC/1.28-foss-2024a

# Directory with reads
reads_dir="/proj/naiss2023-6-65/Milena/chapter3/RNAseq/maculatus/raw_data/fastq"

# Output directory for results
output_dir="/proj/naiss2023-6-65/Milena/chapter3/RNAseq/maculatus/raw_data/trimmed_fastp"
output_reports="/proj/naiss2023-6-65/Milena/chapter3/RNAseq/maculatus/raw_data/trimmed_reports"


for SRR_num in  SRR3113345 SRR3113346 SRR3113354 SRR3113353 SRR3113350 SRR3113352 SRR3113357 SRR3113356 SRR3113361 SRR3113362 SRR3113340 SRR3113341 SRR3113365 SRR3113366 SRR3113378 SRR3113342
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