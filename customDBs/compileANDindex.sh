# Script to pick fasta entries of interest out of the full DB dump of GISAID DB
# The obtained fasta files are merged and outputted as a new fasta file 
# The output fasta file is subsequently used for the indexing of kraken2 DB

#!/bin/bash


targetDir=~/scratch/gisaid
k2FASTAfile=$targetDir/relevantGISAID_k2.fa
kallistoFASTAfile=$targetDir/relevantGISAID_kallisto.fa
taxonomyFile=./majorCovidDB/taxonomy/names.dmp


echo ">wt-wuhan|kraken:taxid|1000" > $k2FASTAfile
cat ../covidRefSequences/wuhan.fa | grep -v '>' >> $k2FASTAfile
echo >> $k2FASTAfile

echo ">wt-wuhan" > $kallistoFASTAfile
cat ../covidRefSequences/wuhan.fa | grep -v '>' >> $kallistoFASTAfile
echo >> $kallistoFASTAfile

echo Generating a subsampled metadata file to improve search time...
gisaidIDs=`cat $taxonomyFile | grep EPI_ISL | awk -F '|' '{print " -e " $2}' | tr -d $'\n' | tr -d $'\t'`
cat $targetDir/metadata.tsv | grep $gisaidIDs > $targetDir/subsampled_metadata.tsv

for gisaidID in `cat $taxonomyFile | grep "EPI_ISL" | awk -F '|' '{print $2}'`; do
	metadataRow=`cat $targetDir/subsampled_metadata.tsv | grep -w $gisaidID`
	fastaHeader=`echo "$metadataRow" | awk -F $'\t' '{ print $1 }'`
	pangoLineage=`echo "$metadataRow" | awk -F $'\t' '{ print $12 }'`
	echo $gisaidID "$fastaHeader" $pangoLineage
	
	grep -A 1000 -m 1 "$fastaHeader" $targetDir/sequences.fasta > $targetDir/excerpt.fa
	head -n 1 $targetDir/excerpt.fa
	csplit -kzq --prefix $targetDir/xx $targetDir/excerpt.fa "/>/" "{*}"
	
	krakenTaxid=`cat $taxonomyFile | grep $gisaidID | awk '{ print $1 }'`
	echo ">$gisaidID|kraken:taxid|$krakenTaxid" >> $k2FASTAfile
	cat $targetDir/xx00 | grep -v '>' >> $k2FASTAfile
	echo >> $k2FASTAfile
	
	echo ">$pangoLineage" >> $kallistoFASTAfile
	cat $targetDir/xx00 | grep -v '>' >> $kallistoFASTAfile
	echo >> $kallistoFASTAfile
	rm $targetDir/xx0*
done


exit 0


# kmer length to use for indexing
#kmer=111
kmer=50

# Index the allCovidDB for kraken2
conda activate c-wap/conda/env-kraken2
#kraken2-build --db allCovidDB --add-to-library $k2FASTAfile
#kraken2-build --build --db allCovidDB --kmer-len $kmer
#bracken-build -d allCovidDB -t 10 -k $kmer

# Index the majorCovidDB for kraken2
../conda/env-kraken2/bin/kraken2-build --db majorCovidDB --add-to-library $k2FASTAfile
../conda/env-kraken2/bin/kraken2-build --build --db majorCovidDB --kmer-len $kmer
../conda/env-kraken2/bin/bracken-build -d majorCovidDB -t 10 -k $kmer
#
conda deactivate


# Generate kallisto index
conda activate c-wap/conda/env-kallisto
../conda/env-kallisto/bin/kallisto index --index variants.kalIdx --make-unique $kallistoFASTAfile 
conda deactivate


