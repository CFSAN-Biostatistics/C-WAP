# Script to pick fasta entries of interest out of the full DB dump of GISAID DB
# Some recent covid reports are selected and aligned against the covid ref sequence.
# The so generated vcf is processed to determine the coordinates of the polymorphic sites.


#!/bin/bash
module load minimap2/2.22 samtools/1.13 ivar


# Select recent samples, with intentional 1% subsampling to reduce workload
gisaid_fa=~/scratch/gisaid/msa_0228/msa_0228.fasta
cat $gisaid_fa | grep -A 1 "00|2022-" | tr -d '-' > subsampled.fa

minimap2 -a --sam-hit-only -2 -x map-hifi ~/c-wap/covidRefSequences/wuhan.mmi subsampled.fa -t 2 -o aligned.sam

# Generation of a sorted bam file from the alignment output
samtools sort aligned.sam -o sorted.bam -@ 2
	
samtools mpileup -aa -A -d 10000 -B -Q 0 --reference ~/c-wap/covidRefSequences/wuhan.fa -o pile.up sorted.bam

varfilename="2022Q1-commonVariants"
cat pile.up | ivar variants -p 2022Q1-commonVariants -g ~/c-wap/covidRefSequences/covidGenomeAnnotation-NCBI.gff \
		-r ~/c-wap/covidRefSequences/wuhan.fa -m 1
		
		
cat ${varfilename}.tsv | grep TRUE | awk '$11 <= 0.1 {print $2}' | uniq > rareSNPpos
cat ${varfilename}.tsv | grep TRUE | awk '$11 > 0.1 && $11 < 0.9 {print $2}' | uniq > diverseSNPpos
cat ${varfilename}.tsv | grep TRUE | awk '$11 >= 0.9 {print $2}' | uniq > commonSNPpos
