#!/bin/bash

# Wastewater covid analysis script. 
# To be executed on a compute node, if available.
# Execution on a local computer is possible, though it might take about an hour.

# Exit if any of these commands fail
set -e


outDir=$1
referenceSequence=$2
primerBedFile=$3
performQConly=$4
platform=$5
sampleName=$6

if [[ -z $8 ]]; then
	singlefilename=$7
else
	R1filename=$7
	R2filename=$8	
fi


echo Analysis workflow started
numThreads=`nproc`
availMem=$(free -g | grep Mem | awk '{ print $2 }') # in GB


# Using half of the available memory just in case, remove it if the RAM available is low.
memPerThread=$(($availMem/(2*$numThreads)))
echo Executing on $numThreads threads with a total memory of $availMem GB, using $memPerThread GB per thread
hostname



#######################################################
###### Alignment and variant calling ##################
#######################################################
echo Aligning reads against the reference...
refSeqBasename=${referenceSequence%.*}

if [[ $platform == Illumina ]]; then
	# Align the reads to the reference sequence to obtain a sorted bam file
	# Check if index files exist, create otherwise
	if ! [[ -f $refSeqBasename.1.bt2 && -f $refSeqBasename.2.bt2 && -f $refSeqBasename.3.bt2 &&
			-f $refSeqBasename.4.bt2 && -f $refSeqBasename.rev.1.bt2 && -f $refSeqBasename.rev.2.bt2 ]]; then
		echo Could not locate Bowtie2 index files for the reference sequence. Rebuilding...
		bowtie2-build $referenceSequence $refSeqBasename
	fi
	
	if [[ -n $singlefilename ]]; then	
		echo Running Bowtie2...
		bowtie2 --no-unal --threads $numThreads -x $refSeqBasename -U $singlefilename -S $outDir/aligned.sam
		
		echo Running kraken2...
		kraken2 $singlefilename --threads $numThreads --report $outDir/k2-std.out > /dev/null
	else	
		echo Running Bowtie2...
		bowtie2 --no-unal --threads $numThreads -x $refSeqBasename -1 $R1filename -2 $R2filename \
			-S $outDir/aligned.sam
		
		echo Running kraken2...
		kraken2 --paired $R1filename $R2filename --threads $numThreads --report $outDir/k2-std.out > /dev/null
	fi
else
	# Align the reads to the reference sequence to obtain a sorted bam file	
	# Check if index file exists, create otherwise
	if ! [[ -f "$refSeqBasename.mmi" ]]; then
		echo Could not locate Minimap2 index files for the reference sequence. Rebuilding...
		minimap2 -d $refSeqBasename.mmi $referenceSequence
	fi
	
	# Perform alignment based on the optimal parameters for the platform
	echo Running minimap2...
	if [[ $platform == Pacbio ]]; then
		minimap2 -a --sam-hit-only -2 -x map-hifi $refSeqBasename.mmi $singlefilename \
			-t $numThreads -o $outDir/aligned.sam
	else
		minimap2 -a --sam-hit-only -2 -x map-ont $refSeqBasename.mmi $singlefilename \
			-t $numThreads -o $outDir/aligned.sam
	fi

	echo Running kraken2...
	kraken2 $singlefilename --threads $numThreads --report $outDir/k2-std.out > /dev/null
fi



#######################################################
# Obtain a pileup file and call variants
echo Sorting and file conversion in progress...
samtools view -bS $outDir/aligned.sam -@ $numThreads > $outDir/aligned.bam
samtools sort $outDir/aligned.bam -o $outDir/sorted.bam -@ $numThreads -m ${memPerThread}G
rm $outDir/aligned.sam $outDir/aligned.bam


# TODO: fix the fragmented/unfragmented library prep issue
echo Starting adaptor trimming...
if [[ $platform == ONT ]]; then
	# Nanopore has a much lower read quality, so the quality trimming should be much more lax.
	ivar trim -b $primerBedFile -p $outDir/trimmed -i $outDir/sorted.bam -q 1
else
	# ivar trim -b $primerBedFile -p $outDir/trimmed -i $outDir/sorted.bam
	# Pass "- e" if fragmented
	cp $outDir/sorted.bam $outDir/trimmed.bam
fi
samtools sort $outDir/trimmed.bam -o $outDir/resorted.bam -@ $numThreads -m ${memPerThread}G
rm $outDir/trimmed.bam


echo Generating a pileup file...
samtools mpileup -aa -A -d 10000 -B -Q 0 --reference $referenceSequence -o $outDir/pile.up $outDir/resorted.bam


echo Calling variants...
minDepthThreshold=10
cat $outDir/pile.up | ivar variants -p $outDir/rawVarCalls -g ./covidRefSequences/covidGenomeAnnotation-NCBI.gff \
	-r $referenceSequence -m $minDepthThreshold



# The pileup file is parsed to calculate the positionwise quality and depth parameters.
# The result is stored as png files that are added to the html report
echo Coverage and quality analysis in progress...
gapfilename=${primerBedFile%.*}.uncovered
./plotCoverageQualityPerPos.py $outDir/pile.up $outDir $gapfilename
# Pile-up file is not needed anymore, so remove to save space.
rm $outDir/pile.up


echo Generating consensus sequence...
# cat $outDir/pile.up | ivar consensus -p $outDir/consensus -t 0.5 -m 1
# Alternative, for cases that the above does not work anymore (segmentation fault):
bcftools mpileup -d 10000 -Ou -f $referenceSequence $outDir/resorted.bam | bcftools call --ploidy 1 -mv -Oz -o $outDir/calls.vcf.gz
bcftools index $outDir/calls.vcf.gz
cat $referenceSequence | bcftools consensus $outDir/calls.vcf.gz > $outDir/consensus.fa
rm $outDir/calls.vcf.gz.csi



#######################################################
###### Variant abundance estimation ###################
#######################################################
if ! $performQConly; then
	#######################################################
	echo Pangolin lineage assignment in progress...
	pangolin --alignment $outDir/consensus.fa --threads $numThreads --outdir $outDir
	mv $outDir/lineage_report.csv $outDir/pangolin_lineage_report.csv 
	rm $outDir/sequences.aln.fasta


	#########################################################################################
	# Generate a new fastq that is covid-only and use for k-mer variant classifiers
	# It is necessary to work with primer-trimmed reads as their presence might deflect the results.
	echo Kraken2 lineage assignment in progress...
	samtools bam2fq $outDir/resorted.bam > $outDir/resorted.fastq

	kraken2 $outDir/resorted.fastq --db ./customDBs/allCovidDB --threads $numThreads\
				--report $outDir/k2-allCovid.out > /dev/null
	bracken -d ./customDBs/allCovidDB -i $outDir/k2-allCovid.out -o $outDir/allCovid.bracken -l P

	kraken2 $outDir/resorted.fastq --db ./customDBs/majorCovidDB --threads $numThreads\
				--report $outDir/k2-majorCovid.out > /dev/null
	bracken -d ./customDBs/majorCovidDB -i $outDir/k2-majorCovid.out -o $outDir/majorCovid.bracken -l P
	rm $outDir/allCovid.bracken $outDir/majorCovid.bracken


	#########################################################################################
	# Generation of index files for all target sequences, if needed:
	# $condaBin/kallisto index --index variants.kalIdx variantGenomeSequences.fa
	echo Running kallisto to quantify the reads hitting variants...
	kallisto quant --index ./customDBs/variants.kalIdx --output-dir $outDir \
			--plaintext -t $numThreads --single -l 500 -s 50 $outDir/resorted.fastq
	rm $outDir/resorted.fastq $outDir/run_info.json
	mv $outDir/abundance.tsv $outDir/kallisto_abundance.tsv
	

	#########################################################################################
	# Variant estimation via Frejya
	echo De-mixing via Freyja...
	freyja variants $outDir/resorted.bam --variants $outDir/freyja.variants.tsv --depths $outDir/freyja.depths --ref $referenceSequence
	freyja demix $outDir/freyja.variants.tsv $outDir/freyja.depths --output $outDir/freyja.demix
	rm $outDir/freyja.variants.tsv $outDir/freyja.depths


	#########################################################################################
	# Linear deconvolution approach
	# Check if the lineage-mutation map file is up-to-date and update otherwise
	variantDBfile=./covidRefSequences/varDefinitions.pkl
	if [ "$variantDBfile" -ot "./constellations/constellations/definitions" ]; then
		echo "WARNING: Pre-processed variant definitions are older than the constellation database. Recompiled."
		mv $variantDBfile "$variantDBfile.old"
		./preprocessVariantDB.py
	else
		echo "Pre-processed variant definitions in $variantDBfile are up-to-date. Using as is."
	fi


	echo Deconvolution-based estimation of the proportion of covid variants
	./deconvolveVariants.py $outDir/rawVarCalls.tsv $outDir $variantDBfile > $outDir/deconvolution.output


	#######################################################
	# Run the relevant python script to generate pie charts figures
	echo Generating abundance plots...
	./plotPieChartsforAbundance.py $outDir $variantDBfile $outDir/linearDeconvolution_abundance.csv \
			$outDir/kallisto_abundance.tsv $outDir/k2-allCovid_bracken.out $outDir/k2-majorCovid_bracken.out $outDir/freyja.demix
fi


echo Done with the analysis on `hostname`.
exit 0

