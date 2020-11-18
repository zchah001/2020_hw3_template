#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --out snp.log

module load samtools
module load bwa
CPU=8


# align the reads from the libraries to the genome
# see bwa mem examples
# libraries are the *.fastq.gz files in your folder now

for libname in $(cat acc.txt)
do
	# YOU NEED TO FIX THIS
	# run the bcftools steps for SNP calling 
done
# after producing the sam file

# caount how many SNPs or INDELs in each file

