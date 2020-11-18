#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --out RNASeq.log


module load kallisto

ln -s /bigdata/gen220/shared/data/M_tuberculosis
ln -s M_tuberculosis/M_tuberculosis.cds.fasta

DB=M_tuberculosis.cds.fasta

# INDEX THE DB FOR KALLISTO
INDEX=$(basename $DB .fasta).idx
# YOUR CODE HERE TO CREATE THIS INDEX

## END INDEX

mkdir -p results


tail -n +2 M_tuberculosis/sra_info.tab | while read RUN CARBON SAMPLE REP
do
	OUTNAME=results/${CARBON}_pH${SAMPLE}_r$REP
	FASTQFILE=M_tuberculosis/$RUN.fastq.gz
	echo "RUN is $RUN $OUTNAME INFILE=$FASTQFILE"
	#RUN KALLISTO ON THE FILE
#	kallisto quant -i $INDEX --single -l 300 -s 30 [ADD SOMETHING HERE]
done
