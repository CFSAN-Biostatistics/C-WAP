#! /usr/bin/env nextflow


PATH="/projects/covidtrakr/software/miniconda3/bin:$PATH"
rootDir="$PWD"

singlefilename = Channel
					.fromPath( "$rootDir/test/test.fastq.gz" )
					.ifEmpty('')
					.into { singlefilename_a; singlefilename_b; singlefilename_c  }

R1filename = Channel
				.empty()
				.ifEmpty('')
				.into { R1filename_a; R1filename_b }

R2filename = Channel
				.empty()
				.ifEmpty('')
				.into { R2filename_a; R2filename_b }


//pairedFilenames = Channel
//    .fromFilePairs('/my/data/SRR*_{1,2}.fastq')
//    .view()


platform = "Illumina"
primerBedFile = "$rootDir/covidRefSequences/none.bed"
variantDBfile = '~/C-WAP/covidRefSequences/varDefinitions.pkl'
referenceSequence = '~/C-WAP/covidRefSequences/wuhan.fa'
performQConly = false
isPairedEnd=false
outdir="/projects/covidtrakr/out"

numThreads = 20


// Align the reads to the reference sequence to obtain a sorted bam file
process runReferenceAlignment {
	// scratch "ram-disk"
	// cpus "$numThreads"
	// memory "20 GB"
	// time "30m"
	
	input:
	file 'single.fastq.gz' from singlefilename_a
	file 'paired_R1.fastq.gz' from R1filename_a
	file 'paired_R2.fastq.gz' from R2filename_a
	
	output:
	file 'resorted.bam' into resorted_bam_a
	file 'resorted.bam' into resorted_bam_b
	file 'resorted.bam' into resorted_bam_c
	file 'resorted.bam' into resorted_bam_d
	file 'sorted.stats' into sorted_stats
	file 'resorted.stats' into resorted_stats

	module 'bowtie2'
	module 'minimap2/2.22'
	module 'samtools/1.13'
	module 'ivar'
	
	"""
	referenceSeqLocal=$referenceSequence
	refSeqBasename=\${referenceSeqLocal%.*}

	case $platform in
		Illumina)
			if [[ -n $singlefilename ]]; then	
				bowtie2 --no-unal --threads $numThreads -x \$refSeqBasename -U single.fastq.gz -S aligned.sam
			else	
				bowtie2 --no-unal --threads $numThreads -x \$refSeqBasename -1 paired_R1.fastq.gz -2 paired_R2.fastq.gz \
					-S aligned.sam			
			fi
			;;
		ONT)
			minimap2 -a --sam-hit-only -2 -x map-ont \${refSeqBasename}.mmi single.fastq.gz \
					-t $numThreads -o aligned.sam
			;;
		Pacbio)
			minimap2 -a --sam-hit-only -2 -x map-hifi \${refSeqBasename}.mmi single.fastq.gz \
					-t $numThreads -o aligned.sam
			;;
	esac
	
	samtools view -bS aligned.sam -@ $numThreads > aligned.bam
	samtools sort aligned.bam -o sorted.bam -@ $numThreads
	
	# Nanopore has a much lower read quality, so the quality trimming should be much more lax.
	if [[ $platform == ONT ]]; then
		ivar trim -e -b $primerBedFile -p trimmed -i sorted.bam -q 1
	else
		ivar trim -e -b $primerBedFile -p trimmed -i sorted.bam
	fi
	
	samtools sort trimmed.bam -o resorted.bam -@ $numThreads
	
	# Evaluate read statistics
	samtools stats sorted.bam | grep ^SN | cut -f 2- > sorted.stats
	samtools stats resorted.bam | grep ^SN | cut -f 2- > resorted.stats
	"""
	// -m ${memPerThread}G
}



process trimmedBam2Fastq {
	input:
	file 'resorted.bam' from resorted_bam_a
	
	output:
	file 'resorted.fastq.gz' into resorted_fastq_gz_a
	file 'resorted.fastq.gz' into resorted_fastq_gz_b
	
	module 'samtools/1.13'
	
	shell:
	"""
	samtools bam2fq resorted.bam > resorted.fastq
	gzip resorted.fastq
	"""
}


process generatePileup {
	input:
	file 'resorted.bam' from resorted_bam_b
	
	output:
	file 'pile.up' into pile_up_a
	file 'pile.up' into pile_up_b


	module 'samtools/1.13'
	
	shell:
	"""
	samtools mpileup -aa -A -d 10000 -B -Q 0 --reference $referenceSequence -o pile.up resorted.bam
	"""
}


process variantCalling {
	input:
	file 'pile.up' from pile_up_a
	
	output:
	file 'rawVarCalls.tsv' into rawVarCalls_tsv
	
	module 'ivar'
	
	shell:
	"""
	cat pile.up | ivar variants -p rawVarCalls -g $rootDir/covidRefSequences/covidGenomeAnnotation-NCBI.gff \
	-r $referenceSequence -m 10
	"""
}


process runKraken2stdDB {
	input:
	file 'singlefilename' from singlefilename_b
	file 'R1filename' from R1filename_b
	file 'R2filename' from R2filename_b
	
	output:
	file 'k2-std.out' into k2_std_out
	
	module '/nfs/software/modules/kraken2/2.1.2'
	
	shell:
	"""	
	if [[ -n $singlefilename ]]; then	
		kraken2 $singlefilename --threads $numThreads --report k2-std.out > /dev/null
	else	
		kraken2 --paired $R1filename $R2filename --threads $numThreads --report k2-std.out > /dev/null		
	fi	
	"""
}


// The pileup file is parsed to calculate the positionwise quality and depth parameters.
// The result is stored as png files that are added to the html report
process plotCoverageQC {
	input:
	file 'pile.up' from pile_up_b
	
	output:
	file 'coverage.png' into coverage_png
	file 'depthHistogram.png' into depthHistogram_png
	file 'quality.png' into quality_png
	file 'qualityHistogram.png' into qualityHistogram_png
	
	file 'pos-coverage-quality.tsv' into pos_coverage_quality_tsv
	
	module 'python/3.8.1'
	
	shell:
	"""
	primerBedFileLocal=$primerBedFile
	gapfilename=\${primerBedFileLocal%.*}.uncovered
	$rootDir/plotCoverageQualityPerPos.py pile.up ./ \$gapfilename
	"""
}




// ////////////////////////////////////////////
// VARIANT CALLING
// ////////////////////////////////////////////

process krakenVariantCaller {
	input:
	file 'resorted.fastq.gz' from resorted_fastq_gz_a
	
	output:
	file 'k2-allCovid_bracken.out' into k2_allCovid_bracken_out_a
	file 'k2-majorCovid_bracken.out' into k2_majorCovid_bracken_out_a
	file 'k2-allCovid.out' into k2_allCovid_out_a
	file 'k2-majorCovid.out' into k2_majorCovid_out_a
	
	file 'k2-allCovid_bracken.out' into k2_allCovid_bracken_out_b
	file 'k2-majorCovid_bracken.out' into k2_majorCovid_bracken_out_b
	file 'k2-allCovid.out' into k2_allCovid_out_b
	file 'k2-majorCovid.out' into k2_majorCovid_out_b
	
	module 'kraken2'
	module 'bracken/2.5'
	
	shell:
	"""
	kraken2 resorted.fastq.gz --db $rootDir/customDBs/allCovidDB --threads $numThreads --report k2-allCovid.out > /dev/null
	bracken -d $rootDir/customDBs/allCovidDB -i k2-allCovid.out -o allCovid.bracken -l P

	kraken2 resorted.fastq.gz --db $rootDir/customDBs/majorCovidDB --threads $numThreads --report k2-majorCovid.out > /dev/null
	bracken -d $rootDir/customDBs/majorCovidDB -i k2-majorCovid.out -o majorCovid.bracken -l P
	"""
}



process pangolinVariantCaller {
	input:
	file 'resorted.bam' from resorted_bam_c
	
	output:
	file 'lineage_report.csv' into pangolin_lineage_report_csv
	val 'consensusLineage' into consensusLineage
	
	module 'samtools/1.13'
	conda '/projects/covidtrakr/software/miniconda3/envs/pangolin'
	
	shell:
	"""	
	bcftools mpileup -d 10000 -Ou -f $referenceSequence resorted.bam | bcftools call --ploidy 1 -mv -Oz -o calls.vcf.gz
	bcftools index calls.vcf.gz
	cat $referenceSequence | bcftools consensus calls.vcf.gz > consensus.fa
	
	pangolin --alignment consensus.fa --threads $numThreads --outdir ./
	
	# Characterisation of the consensus sequence based on Pangolin output
	# Calculation of the consensus sequence is used to determine the predominant lineage.
	# If available, use the WHO label.
	consensusLineage=`tail -n 1 pangolin_lineage_report.csv | awk -F "," '{ print \$2 }'`
	WHOlabel=`cat $rootDir/pangolin2WHOlabel.txt | grep \$consensusLineage | awk -F " " '{ print \$2 }'`
	if  [[ -n \$WHOlabel ]]; then
		consensusLineage=\$WHOlabel
	fi
	"""
}



process linearDeconVariantCaller {
	input:
	file 'rawVarCalls.tsv' from rawVarCalls_tsv
	
	output:
	file 'linearDeconvolution_abundance.csv' into linearDeconvolution_abundance_csv_a
	file 'linearDeconvolution_abundance.csv' into linearDeconvolution_abundance_csv_b
	file 'mutationTable.html' into mutationTable_html
	file 'VOC-VOIsupportTable.html' into VOC_VOIsupportTable_html
	
	val 'mostAbundantVariantPct' into linRegressionTopPct
	val 'mostAbundantVariantName' into linRegressionTopName
	val 'linRegressionR2' into linRegressionR2
	
	module 'python/3.8.1'
	
	shell:
	"""
	deconvolutionOutput=`$rootDir/deconvolveVariants.py rawVarCalls.tsv ./ $variantDBfile`
	
	mostAbundantVariantPct=`echo \$deconvolutionOutput | awk '{ print \$1 }'`
	mostAbundantVariantName=`echo \$deconvolutionOutput |awk '{ print \$2 }'`
	linRegressionR2=`echo \$deconvolutionOutput | awk '{ print \$3 }'`
	"""
}


process kallistoVariantCaller {
	input:
	file 'resorted.fastq.gz' from resorted_fastq_gz_b
		
	output:
	file 'abundance.tsv' into kallisto_abundance_tsv_a
	file 'abundance.tsv' into kallisto_abundance_tsv_b
	
	conda 'kallisto'
	
	shell:
	"""
	numThreads=20
	kallisto quant --index $rootDir/customDBs/variants.kalIdx --output-dir ./ \
			--plaintext -t $numThreads --single -l 500 -s 50 resorted.fastq.gz	
	"""
}


process freyjaVariantCaller {
	input:
	file 'resorted.bam' from resorted_bam_d
		
	output:
	file 'freyja.demix' into freyja_demix_a
	file 'freyja.demix' into freyja_demix_b

	conda '/projects/covidtrakr/software/miniconda3/envs/freyja-env'
	
	// Due to a potential bug, some big fastqs result in a pandas error.
	// Start by generating an empty file to circumvent such failure cases
	
	shell:
	"""
	echo FATAL ERROR > freyja.demix
	echo summarized\$'\t'"[('Other', 1.00)]" >> freyja.demix
	
	echo Pileup generation for Freyja
	freyja variants resorted.bam --variants freyja.variants.tsv --depths freyja.depths --ref $referenceSequence || true
	
	echo Demixing variants by Freyja
	freyja demix freyja.variants.tsv freyja.depths --output freyja.demix || true
	"""
}



///////////////////////////////////////////////
// Report generation and final output
///////////////////////////////////////////////

process plotPieCharts {
	input:
	file 'linearDeconvolution_abundance.csv' from linearDeconvolution_abundance_csv_a
	file 'kallisto_abundance.tsv' from kallisto_abundance_tsv_a
	file 'k2-allCovid_bracken.out' from k2_allCovid_bracken_out_a
	file 'k2-majorCovid_bracken.out' from k2_majorCovid_bracken_out_a
	file 'freyja.demix' from freyja_demix_a
	
	output:
	file 'pieChart_k2_allCovid.png' into pieChart_k2_allCovid_png
	file 'pieChart_k2_majorCovid.png' into pieChart_k2_majorCovid_png
	file 'pieChart_kallisto.png' into pieChart_kallisto_png
	file 'pieChart_deconvolution.png' into pieChart_deconvolution_png
	file 'pieChart_freyja.png' into pieChart_freyja_png

	val 'kallistoTopName' into kallistoTopName

	module 'python/3.8.1'
	
	shell:
	"""
	$rootDir/plotPieChartsforAbundance.py ./ $variantDBfile linearDeconvolution_abundance.csv \
			kallisto_abundance.tsv k2-allCovid_bracken.out k2-majorCovid_bracken.out freyja.demix

	kallistoTopName=`cat kallisto.out | sort -k 2 -n | tail -n 1 | awk '{ print \$1 }'`
	"""
}



// A metadata fetch attemp from NCBI via Entrez-Direct
// Only works if file names explicitly carries an SRR number.
process getNCBImetadata {
	input:
	file 'singlefilename' from singlefilename_c
	
	output:
	val libraryProtocol
	val seqInstrument
	val isolate
	val collectionDate
	val collectedBy
	val sequencedBy
	val sampleLatitude
	val sampleLongitude
	val sampleLocation
	
	module 'edirect'
	
	shell:
	"""
	SAMN=Missing
	libraryProtocol=Missing
	seqInstrument=Missing
	isolate=Missing
	collectionDate=Missing
	collectedBy=Missing
	sequencedBy=Missing
	sampleLatitude=0
	sampleLongitude=0
	sampleLocation=NA
	
	sfilename=$singlefilename
	srrNumber=\${sfilename%_*}
	if [[ \${srrNumber:0:3} == 'SRR' ]]; then
		# The tool returns error: too many requests, bypassing by redirection of error
		sraQueryResult=\$(esearch -db sra -query \$srrNumber 2>/dev/null)
		if echo \$sraQueryResult | grep -q "<Count>1</Count>"; then
			# Get runinfo from SRA
			echo Downloading and parsing metadata for \$srrNumber...	
			echo "\$sraQueryResult" | efetch --format runinfo
			SRRmetadata=`echo "\$sraQueryResult" | efetch --format runinfo 2>/dev/null | grep \$srrNumber`
			echo 000 \$SRRmetadata
			
			libraryProtocol=`echo \$SRRmetadata | awk -F ',' '{print \$13}'`
			seqInstrument=`echo \$SRRmetadata | awk -F ',' '{print \$20}'`
			isolate=`echo \$SRRmetadata  | awk -F ',' '{print \$30}'`	
			
			# Get metadata out of biosample db
			echo Fetching biosample data
			SAMN=`echo \$SRRmetadata | awk -F ',' '{print \$26}'`
			SAMNmetadata=`efetch -db biosample -id \$SAMN 2>/dev/null`
			
			collectionDate=`echo "\$SAMNmetadata" | grep "collection date" | awk -F '"' '{print \$2}'`
			collectedBy=`echo "\$SAMNmetadata" | grep "collected by" | awk -F '"' '{print \$2}'`
			sequencedBy=`echo "\$SAMNmetadata" | grep SEQUENCED_BY | awk '{ \$1=""; print \$0 }'`
			sampleLatitude=`echo "\$SAMNmetadata" | grep "latitude and longitude" | awk -F '"' '{print \$2}'\
								| awk '{ print \$1\$2 }'`
			sampleLongitude=`echo "\$SAMNmetadata" | grep "latitude and longitude" | awk -F '"' '{print \$2}'\
								| awk '{ print \$3\$4 }'`
			sampleLocation=`echo "\$SAMNmetadata" | grep "geographic location" | awk -F '"' '{print \$2}'`
		fi
	fi

	if [[ -z "\$SAMN" ]]; then
		echo No SRA metadata was available for \$srrNumber
	fi
	"""
}



process generateReport {
	publishDir "$outdir", mode: 'copy', overwrite: true
	
	input:
	file 'k2-std.out' from k2_std_out
	file 'mutationTable.html' from mutationTable_html
	file 'VOC-VOIsupportTable.html' from VOC_VOIsupportTable_html
	
	val libraryProtocol
	val seqInstrument
	val isolate
	val collectionDate
	val collectedBy
	val sequencedBy
	val sampleLatitude
	val sampleLongitude
	val sampleLocation
	
	file 'pieChart_k2_allCovid.png' from pieChart_k2_allCovid_png
	file 'pieChart_k2_majorCovid.png' from pieChart_k2_majorCovid_png
	file 'pieChart_kallisto.png' from pieChart_kallisto_png
	file 'pieChart_deconvolution.png' from pieChart_deconvolution_png
	file 'pieChart_freyja.png' from pieChart_freyja_png
	
	file 'sorted.stats' from sorted_stats
	file 'resorted.stats' from resorted_stats

	file 'pos-coverage-quality.tsv' from pos_coverage_quality_tsv
	file 'pangolin_lineage_report.csv' from pangolin_lineage_report_csv
	
	env 'mostAbundantVariantPct' from linRegressionTopPct
	env 'mostAbundantVariantName' from linRegressionTopName
	env 'linRegressionR2' from linRegressionR2
	env 'kallistoTopName' from kallistoTopName
	
	file 'coverage.png' from coverage_png
	file 'depthHistogram.png' from depthHistogram_png
	file 'quality.png' from quality_png
	file 'qualityHistogram.png' from qualityHistogram_png
	
	file 'linearDeconvolution_abundance.csv' from linearDeconvolution_abundance_csv_b
	file 'kallisto_abundance.tsv' from kallisto_abundance_tsv_b
	file 'k2-allCovid_bracken.out' from k2_allCovid_bracken_out_b
	file 'k2-majorCovid_bracken.out' from k2_majorCovid_bracken_out_b
	file 'freyja.demix' from freyja_demix_b

	output:
	file 'a'
	// file '*.stats'
	//file '*.csv'
	//file '*.tsv'
	file 'k2*n.out' into k2Out
	file 'freyja.demix'
	
	shell:
	"""
	cp $rootDir/htmlHeader.html ./
	
	sfilename=$singlefilename
	sampleName=\${sfilename%_*}
	$rootDir/generateReport_nf.sh \$sampleName $primerBedFile $performQConly $isPairedEnd

	wkhtmltopdf --enable-local-file-access --page-size Letter --margin-top 10mm --margin-bottom 0 \
			--margin-left 0 --margin-right 0 --print-media-type --title "Wastewater report" \
			report.html report.pdf
	
	mkdir report
	mv report.* report/
	mv *.png report/
	
	mkdir a
	mv * a/
	"""
}



process summaryPage {
	input:

	output:
	
	module 'bowtie2'
	module 'minimap2'
	module 'samtools/1.13'
	conda 'kraken2'
	conda 'kallisto'
	conda '/projects/covidtrakr/software/miniconda3/envs/freyja-env'
	conda '/projects/covidtrakr/software/miniconda3/envs/pangolin'


	shell:
	"""
	reportFile="report.html"
	echo >> \$reportFile
	
	# Appending a list of data analysis parameters...
	echo "<h2>Software configuration</h2>" >> \$reportFile

	pangolinVersion=`pangolin -v`
	pangolearnVersion=`pangolin -pv`

	if [[ -n "$singlefilename" ]]; then
		echo "Minimap2:" `minimap2 --version` >> \$reportFile
	else
		echo "Bowtie2 v`bowtie2 --version | head -n 1 | awk '{ print \$3 }'`, " >> \$reportFile
	fi

	echo `samtools --version | head -n 2`", " >> \$reportFile
	echo `ivar version | head -n 1`", " >> \$reportFile
	echo `kraken2 -v | head -n 1`", " >> \$reportFile
	echo `kallisto version`"." >> \$reportFile

	allIncludedLineages=`$rootDir/listVariantsAvail.py $variantDBfile`
	echo "Lineage definitions were compiled on" `date +%Y-%m-%d -r ./covidRefSequences/varDefinitions.pkl` \
			"from <a href=\"https://github.com/cov-lineages/constellations/tree/main/constellations/definitions\">constellations</a>." >> \$reportFile
	echo "Lineage signature file was compiled on" `date +%Y-%m-%d -r \$variantDBfile` \
			"and includes lineages: \$allIncludedLineages." >> \$reportFile
	echo "Lineage assignment to the consensus sequence was performed by \$pangolinVersion using the classification tree of \$pangolearnVersion." >> \$reportFile
	"""
}

