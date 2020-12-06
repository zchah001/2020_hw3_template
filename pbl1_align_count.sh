#!/usr/bin/bash
#SBATCH -p short -N 1 -n 4 --mem 8gb --out make_links.log

module load bwa
module load samtools

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
ACC="SRR11587604_1 SRR11587604_2 SRR11140748_1 SRR11140748_2"
for acc in $ACC;
do
  	echo "$acc"
        if [ ! -s ${acc}_1.fastq.gz ]; then
                ln -s $FOLDER/${acc}_[ZC].fastq.gz .
        fi
done


#!/usr/bin/bash
#SBATCH -p batch -N 1 -n 16 --mem 8gb --out align.log

module load bwa
module load samtools
CPU=16

# make a link
GENOME=NC_045512.fa

ACC="SRR11587604_1 SRR11587604_2 SRR11140748_1 SRR11140748_2"
for acc in $ACC;
do
  	if [ ! -s ${acc}.fixmate.bam ]; then
                bwa mem -t $CPU $GENOME ${acc}_[ZC].fastq.gz > $acc.sam
                samtools fixmate --threads $CPU -O bam $acc.sam ${acc}.fixmate.bam
        fi
	if [ ! -s ${acc}.bam ]; then
                samtools sort -O bam --threads $CPU -o ${acc}.bam ${acc}.fixmate.bam
        fi
	if [ -f ${acc}.bam ]; then
                samtools flagstat ${acc}.bam > ${acc}.stats.txt
        fi
done


#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 8gb --out align.%a.log

module load bwa
module load samtools
CPU=8
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi

# make a link
GENOME=NC_045512.fa
ACCFILE=acc.txt
sed -n ${N}p $ACCFILE | while read acc;
do
  	if [ ! -f ${acc}.arrayjobs.fixmate.bam ]; then
                bwa mem -t $CPU $GENOME ${acc}[Z].fastq.gz  > ${acc}.arrayjobs.sam
                samtools fixmate --threads $CPU -O bam ${acc}.arrayjobs.sam ${acc}.arrayjobs.fixmate.bam
        fi
	if [ ! -f ${acc}.arrayjobs.bam ]; then
                samtools sort -O bam --threads $CPU -o ${acc}.arrayjobs.bam ${acc}.arrayjobs.fixmate.bam
        fi
	if [ ! -f ${acc}.arrayjobs.stats.txt ]; then
                samtools flagstat ${acc}.arrayjobs.bam > ${acc}.arrayjobs.stats.txt
        fi
done
