#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -c 12 
#SBATCH --mem=100G
#SBATCH -t 5:00:00
#SBATCH -J SNP_calling
#SBATCH -o SNP_calling.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se


### This is from Jean-Loup
### run like this:
#sbatch 3_GATK.sh -b /path/to/bam/bamfile.bam -r /path/to/bedfile/bedfile -i samplename -o outdir -g /path/to/refgenome/refgenome.fa


module unload
module load SAMtools/1.22-GCC-13.3.0
module load GATK/4.6.2.0-GCCcore-13.3.0-Java-17


while getopts b:r:i:o:g: flag 
do
case "${flag}" in
        b) BAM=${OPTARG}
            ;;
        o) outdir=${OPTARG}
            ;;
        g) ref=${OPTARG}
            ;;
        *) echo "Invalid option: -$flag" 
            ;;
    esac
done

echo "## ------------------------------------------------ ##"
echo `date` ": sync files."
rsync -ta $BAM $TMPDIR
rsync -ta $ref* $TMPDIR
dict="${ref_prefix}.dict"
echo "## ------------------------------------------------ ##"
echo "Initiating program"

cd $TMPDIR
echo
ls $TMPDIR
echo

genome=`echo $ref |awk -F '/' '{print $NF}'`
bam=`echo $BAM |awk -F '/' '{print $NF}'`

echo "looking at $bam"
echo "Reference is $genome" 

if [[ ! -f "${bam}.bai" ]]; then
    echo "Index missing in TMP, indexing BAM"
    samtools index "$bam"
fi

dict=${ref}.dict
if [[ ! -f "$dict" ]]; then
    gatk --java-options "-Xmx32G" CreateSequenceDictionary -R $genome
fi

echo
ls $TMPDIR
gatk --java-options "-Xmx91G" HaplotypeCaller -R $genome -I $bam -ERC GVCF --tmp-dir $TMPDIR -O ${bam%.bam}.vcf.gz 


mkdir -p $outdir
rsync -r ${bam%.bam}.vcf.gz $outdir
