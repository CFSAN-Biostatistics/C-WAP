#!/bin/bash

# Script to pick fasta entries of interest out of the full DB dump of GISAID DB

# Download a copy of the GISAID sequences (full unmasked) before starting, and extract.


# Some recent covid reports are selected and aligned against the covid ref sequence.
# The so generated vcf is processed to determine the coordinates of the polymorphic sites.


module load minimap2/2.22 samtools/1.13 ivar


#GISAIDdir=~/scratch/GISAID/
GISAIDdir=/hpc/scratch/$USER/GISAID/


# Select recent samples, with intentional 1% subsampling to reduce workload
# Subsampling is done by only retaining submissions whose IDs end with "00".
echo Subsampling...
#gisaid_fa=~/scratch/GISAID/msa_0721.fasta
gisaid_fa=/hpc/scratch/$USER/GISAID/msa_0809/msa_0809.fasta
cat $gisaid_fa | grep -A 1 "00|2022-" | tr -d '-' > $GISAIDdir/subsampled.fa

echo Aligning...
minimap2 -a --sam-hit-only -2 -x map-hifi ../covidRefSequences/wuhan.mmi $GISAIDdir/subsampled.fa -t 20 -o $GISAIDdir/aligned.sam

# Generation of a sorted bam file from the alignment output
echo Sorting...
samtools sort $GISAIDdir/aligned.sam -o $GISAIDdir/sorted.bam -@ 20

echo Pileup generation in progress...
samtools mpileup -aa -A -d 10000 -B -Q 0 --reference ../covidRefSequences/wuhan.fa -o $GISAIDdir/pile.up $GISAIDdir/sorted.bam

echo Variant calls...
varfilename=commonVariants
cat $GISAIDdir/pile.up | ivar variants -p $varfilename -g ../covidRefSequences/covidGenomeAnnotation-NCBI.gff \
		-r ../covidRefSequences/wuhan.fa -m 1

cat ${varfilename}.tsv | grep TRUE | awk '$11 <= 0.1 {print $2}' | uniq > rareSNPpos
cat ${varfilename}.tsv | grep TRUE | awk '$11 > 0.1 && $11 < 0.9 {print $2}' | uniq > diverseSNPpos
cat ${varfilename}.tsv | grep TRUE | awk '$11 >= 0.9 {print $2}' | uniq > commonSNPpos

echo Done.
