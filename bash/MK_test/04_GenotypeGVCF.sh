#!/bin/bash
#SBATCH -J genotypeGVCF
#SBATCH -A uppmax2026-1-8
#SBATCH -c 16 
#SBATCH -t 2-23:59:00

### This is from Jean-Loup
### run like this:
## sbatch GenotypeGVCF.sh -d /path/to/database/ -o /path/to/outdir -n outfilename -g path/to/refgenome 

ulimit -c unlimited

module load GATK/4.6.2.0-GCCcore-13.3.0-Java-17
module load BCFtools/1.22-GCC-13.3.0


while getopts d:o:n:g: flag
do
case "${flag}" in
                d) DIRDATA=${OPTARG}
                    	;;
                o) outdir=${OPTARG}
                        ;;
                n) name=${OPTARG}
			            ;;
                g) ref=${OPTARG}
			            ;;
                *) echo "Invalid option: -$flag" 
                      	;;


        esac
done


echo "## ------------------------------------------------ ##"
echo `date` ": sync files."
cp -r $DIRDATA $TMPDIR
cp $ref* $TMPDIR
cd $TMPDIR

genome=`echo $ref |awk -F '/' '{print $NF}'`

if [[ ! -f "${genome}.dict" ]]; then
    gatk --java-options "-Djava.io.tmpdir=tmp -Xmx32G"  CreateSequenceDictionary -R $genome -O ${genome%.fa}.dict
fi


gatk --java-options "-Xmx110G" GenotypeGVCFs -V gendb://my_database -R $genome -O $name.vcf


bgzip "$name".vcf
cd -

mkdir -p $outdir
cp $TMPDIR/$name* $outdir
