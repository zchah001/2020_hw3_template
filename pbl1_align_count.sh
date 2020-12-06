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

# make a link
FOLDER=/bigdata/gen220/shared/data/SARS-CoV-2/
GENOME=NC_045512.fa
if [ ! -L $GENOME ]; then
        # make a symlink
        ln -s $FOLDER/$GENOME .
fi

if [ ! -f $GENOME.pac ]; then
        bwa index $GENOME
fi
acc="SRR11587604_1 SRR11587604_2 SRR11140748_1 SRR11140748_2"
for acc in $acc;
do
        echo "$acc"
        if [ ! -s ${acc}_1.fastq.gz ]; then
                ln -s $FOLDER/${acc}_[Z12].fastq.gz .
        fi



# align the reads from the libraries to the genome
# see bwa mem examples
# libraries are the *.fastq.gz files in your folder now

for acc in $(cat acc.txt)
    do
    # YOU NEED TO FIX THIS
    if [ ! -f $acc.sam ]; then
        bwa mem -t $CPU -o ${acc}.sam $GENOME ${acc$}
    fi
    samtools fixmate --threads $CPU $acc.sam $acc$
    # libraries are the *.fastq.gz files in your folder now

for acc in $(cat acc.txt)
    do
    # YOU NEED TO FIX THIS
    if [ ! -f $acc.sam ]; then
        bwa mem -t $CPU -o ${acc}.sam $GENOME ${acc$}
    fi
    samtools fixmate --threads $CPU $acc.sam $acc$
    samtools sort --threads $CPU -o $acc.bam $acc$
    samtools flagstat $acc.bam > $acc$.stats.txt
done

