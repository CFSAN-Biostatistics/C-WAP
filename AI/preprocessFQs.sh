#! /bin/bash

source ~/c-wap/CFSANonly/prepareEnvironment.sh


numThreads=15

for fqfile in $(ls fastq/*_R1*.gz); do
	echo Now processing $fqfile...
	bowtie2 --no-unal --threads $numThreads -x ~/c-wap/covidRefSequences/wuhan -1 $fqfile -2 $(echo $fqfile | sed 's/\_R1/\_R2/') -S aligned.sam
	samtools sort aligned.sam -o sorted.bam -@ $numThreads

	freyja variants sorted.bam --variants freyja.variants --depths freyja.depths --ref ~/c-wap/covidRefSequences/wuhan.fa
	freyja demix freyja.variants.tsv freyja.depths --output freyja.demix

	rm aligned.sam sorted.bam
	fqname=$(basename $fqfile | awk -F '_' '{print $1}')
	mkdir inputs/$fqname
	mv freyja.demix inputs/$fqname/
	mv freyja.variants.tsv inputs/$fqname/freyja.variants
	mv freyja.depths inputs/$fqname/
done

