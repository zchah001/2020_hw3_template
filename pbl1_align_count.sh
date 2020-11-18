#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --out bwa_samtools.log

module load samtools
module load bwa
CPU=8
ln -s /bigdata/gen220/shared/data/SARS-CoV-2
ln -s SARS-CoV-2/NC_045512.fa .
ln -s SARS-CoV-2/*.fastq.gz .
ln -s SARS-CoV-2/acc.txt

# INDEX the GENOME with `bwa index`

# YOUR CODE 


# align the reads from the libraries to the genome
# see bwa mem examples
# libraries are the *.fastq.gz files in your folder now

for libname in $(cat acc.txt)
do
	# YOU NEED TO FIX THIS
	bwa mem -t $CPU -o [OUTNAME] [GENOME FILE] [FORWARD READ] [REVERSE READ]
done
# after producing the sam file

samtools view -O bam -o [outbame] [insam]

samtools sort [fixme]

samtools flagstat [fixme]
