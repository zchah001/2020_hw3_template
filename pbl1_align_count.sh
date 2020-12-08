#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --out bwa_samtools.log

module load samtools
module load bwa
CPU=8
ln -s /bigdata/gen220/shared/data/SARS-CoV-2
ln -s SARS-CoV-2/NC_045512.fa .
ln -s SARS-CoV-2/*.fastq.gz .
ln -s SARS-CoV-2/acc.txt

GENOME=NC_045512.fa
# INDEX the GENOME with `bwa index`
if [ ! -f $GENOME.pac ]; then
    bwa index $GENOME
fi
# YOUR CODE
libname="SRR11587604 SRR11140742"

# align the reads from the libraries to the genome
# see bwa mem examples
# libraries are the *.fastq.gz files in your folder now

for libname in $(cat acc.txt)
do
    # YOU NEED TO FIX THIS
    if [ ! -f $libname.sam ]; then
        bwa mem -t $CPU -o ${libname}.sam $GENOME ${libname}_1.fastq.gz ${libname}_2.fastq.gz
    fi
    samtools fixmate --threads $CPU $libname.sam $libname.fixmate.bam
    samtools sort --threads $CPU -o $libname.bam $libname.fixmate.bam
    samtools flagstat $libname.bam > $libname.stats.txt
done
