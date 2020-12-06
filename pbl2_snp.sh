#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --out snp.log

module load samtools
module load bwa
module load bcftools
CPU=8

FOLDER=/bigdata/gen220/shared/data/SARS-CoV-2/
GENOME=NC_045512.fa
bwa index $GENOME

for acc  in $(cat acc.txt)
do
  	m="${acc}.bam $m"
        
VCF=SARS-CoV-2.vcf.gz
bcftools mpileup -Ou -f $GENOME $m | bcftools call --ploidy 1 -vmO z -o $VCF
tabix -p vcf $VCF
bcftools stats -F $GENOME -s - $VCF > $VCF.stats
mkdir -p plots
plot-vcfstats -p plots/ $VCF.stats

done
