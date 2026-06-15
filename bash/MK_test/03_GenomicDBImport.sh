#!/bin/bash
#SBATCH -J JoinCalling_add
#SBATCH -A uppmax2026-1-8
#SBATCH -c 16
#SBATCH -t 23:59:00
##cree la database avant le JoinCalling des gvcf 

### This is from Jean-Loup
### run like this:
##sbatch GenomicDBImport.sh -d path/to/dir/with/singlesample.vcf -o path/to/outdir -r path/to/referencegenome 

set -euo pipefail
ulimit -c unlimited

module load SAMtools/1.22-GCC-13.3.0
module load GATK/4.6.2.0-GCCcore-13.3.0-Java-17

while getopts d:o:r: flag; do
    case "${flag}" in
        d) DIRDATA=${OPTARG} ;;
        o) outdir=${OPTARG} ;;
        r) ref=${OPTARG} ;;
        *) echo "Invalid option: -$flag"; exit 1 ;;
    esac
done

echo "## ------------------------------------------------ ##"
echo "$(date): syncing files."
rsync -ta "$DIRDATA"/* "$TMPDIR"/
rsync -ta "${ref}"* "$TMPDIR"/

cd "$TMPDIR"
genome=$(basename "$ref")

for i in *.vcf.gz; do
    [[ ! -f "${i}.tbi" ]] && gatk --java-options "-Xmx80G" IndexFeatureFile -I "$i"
    printf '%s\t%s\n' "${i%.vcf.gz}" "./${i}" >> cohort.sample_map
done

# Refindexing
[[ ! -f "${genome}.fai" ]]  && samtools faidx "$genome"
[[ ! -f "${genome%.fa*}.dict" ]] && gatk CreateSequenceDictionary -R "$genome"

awk '{print $1"\t1\t"$2}' "${genome}.fai" > intervals.list


gatk --java-options "-Xmx100G" GenomicsDBImport \
    --genomicsdb-workspace-path ./my_database \
    --batch-size 20 \
    --sample-name-map cohort.sample_map \
    -L intervals.list

cd -

mkdir -p "$outdir"
rsync -r "$TMPDIR"/my_database "$outdir"/
rsync "$TMPDIR"/*.tbi "$outdir"/ 2>/dev/null || echo "Warning: no .tbi files to copy"