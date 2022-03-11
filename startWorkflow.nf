#! /usr/bin/env nextflow


switch(params.platform) {
	case "i":
	case "Illumina":
		print("Paired end Illumina mode")
		platform = "Illumina"
		isPairedEnd = true
		break;
	case "s":
		print("Single read Illumina mode")
		platform = "Illumina"
		isPairedEnd = false
		break;
	case "n":
	case "ONT":
	case "nanopore":
		print("Nanopore mode")
		platform = "ONT"
		isPairedEnd = false
		break;
	case "p":
	case "pb":
	case "PB":
	case "pacbio":
	case "PacBio":
		print("PacBio mode")
		platform = "PacBio"
		isPairedEnd = false
		break;
	case "h":
		printUsage()
		return 0
	default:
		print("ERROR: unknown sequencing modality")
		printUsage()
		return 1
}


// Convert relative path to absolute path
if (params.primerBedFile[0]=='/' || params.primerBedFile[0]=='~')
	primerBedFile=params.primerBedFile
else
	primerBedFile="$launchDir/$params.primerBedFile"


def printUsage () {
	println("Example usage: ")
	println("\t./startWorkflow.nf --platform i --primers path/to/bed --in path/to/fastq/ --out path/to/outputDir")
	println()
	println("Run the C-WAP workflow that will import the fastq files, trim adaptors, apply quality trimming, call variants and generate a html and pdf report and quit. All fastq files in the provided directory will be processed.")
	println()
	
	println("Platform options:")
	println("\t-h\t\t: Show this help message and exit.")
	println("\t-i\t\t: Use parameters optimised for the Illumina platform, paired-end mode")
	println("\t-s\t\t: Use parameters optimised for the Illumina platform, single read mode")
	println("\t-n\t\t: Use parameters optimised for the ONT platform")
	println("\t-p\t\t: Use parameters optimised for PacBio platform")
}


// Given a file name, extracts a human readable sample name to be used in the output.
allSampleNames = []
def getSampleName(filename) {
	// Typical for most users using Illumina output as is
	// Ex: /path/to/dir/something_S1_L2_R1.fastq -> something
	sampleName = filename.name.split("/")[-1].split("\\.")[0].split("_")[0]
	
	// For complicated file names involving underscores that cannot be eliminated
	// Ex: /path/to/dir/some_thing_R1.fastq -> some-thing	
	// sampleName = filename.name.split("/")[-1].split("\\.")[0].split("_R")[0].replace('_','-')
	
	// If there are special fixed substrings available within all file names.
	// Ex: /path/to/dir/some_thing_1art_out_R1.fastq -> some-thing-1
	// sampleName = filename.name.split("/")[-1].split("\\.")[0].split("art_out")[0].replace('_','-')
	
	// Check if there is name collision and immediately abort the execution if this is the case.
	if (allSampleNames.contains(sampleName)) {
		throw new IllegalArgumentException ('The sample name is not unique, please adjust the getSampleName function.')
	}
	else {
		allSampleNames += sampleName
		return sampleName
	}
}


// Import the list of files to process
if (isPairedEnd) {
	FQs = Channel
	    .fromFilePairs("$params.in/*_R{1,2}*.fastq*", checkIfExists: true, flat:true)
		.map{ tuple(getSampleName(it[1]), it[1], it[2]) }
}
else {
	FQs = Channel
		.fromPath( "$params.in/*.fastq*", checkIfExists: true)
		.map{ tuple(getSampleName(it), it, null) }
}


FQs
	.view()
	.into{ input_fq_a; input_fq_b; input_fq_c; input_fq_d; input_fq_e }



//////////////////////////////////////////////////////////////////////////
// Read processing using common tools, QC etc.
//////////////////////////////////////////////////////////////////////////

// Align the reads to the reference sequence to obtain a sorted bam file
process referenceAlignment {
	cpus 10
	memory '16 GB'
	time = '1 h'

	input:
		tuple val(sampleName), file('R1.fastq.gz'), file('R2.fastq.gz') from input_fq_a
	
	output:
		tuple val(sampleName), env(numReads), path('resorted.bam') into resorted_bam_a, resorted_bam_b,  resorted_bam_c,  resorted_bam_d
		tuple val(sampleName), path('sorted.stats'), path('resorted.stats') into samtools_stats
	
	shell:
	refSeqBasename = params.referenceSequence.replaceAll('.fa$', '')
	"""
	numThreads=\$(nproc)
	case $platform in
		Illumina)
			if $isPairedEnd; then
				bowtie2 --no-unal --threads \$numThreads -x $refSeqBasename -1 R1.fastq.gz -2 R2.fastq.gz \
					-S aligned.sam
			else
				bowtie2 --no-unal --threads \$numThreads -x $refSeqBasename -U R1.fastq.gz -S aligned.sam
			fi
			;;
		ONT)
			minimap2 -a --sam-hit-only -2 -x map-ont ${refSeqBasename}.mmi R1.fastq.gz \
					-t \$numThreads -o aligned.sam
			;;
		PacBio)
			minimap2 -a --sam-hit-only -2 -x map-hifi ${refSeqBasename}.mmi R1.fastq.gz \
					-t \$numThreads -o aligned.sam
			;;
	esac
	
	# Generation of a sorted bam file from the alignment output
	samtools sort aligned.sam -o sorted.bam -@ \$numThreads
	rm aligned.sam
	
	# Nanopore has a much lower read quality, so the quality trimming should be much more lax.
	if [[ $platform == ONT ]]; then
		ivar trim -e -b $primerBedFile -p trimmed -i sorted.bam -q 1
	else
		ivar trim -e -b $primerBedFile -p trimmed -i sorted.bam
	fi
	
	samtools sort trimmed.bam -o resorted.bam -@ \$numThreads
	
	# Evaluate read statistics
	samtools stats sorted.bam | grep ^SN | cut -f 2- > sorted.stats
	samtools stats resorted.bam | grep ^SN | cut -f 2- > resorted.stats
	numReads=\$(cat resorted.stats | grep "raw total sequences" | awk '{ print \$4 }')
	"""
}


process trimmedBam2Fastq {
	input:
		tuple val(sampleName), env(numReads), path('resorted.bam') from resorted_bam_a
	
	output:
		tuple val(sampleName), env(numReads), path('resorted.fastq.gz') into resorted_fastq_gz_a, resorted_fastq_gz_b, resorted_fastq_gz_c
		
	shell:
	"""
		# include -h for header
		samtools view --threads 2 resorted.bam -o resorted.sam
		$projectDir/sam2fastq.py resorted.sam resorted.fastq
		gzip resorted.fastq
		rm resorted.sam
	"""
}


process generatePileup {
	input:
		tuple val(sampleName), env(numReads), path('resorted.bam') from resorted_bam_b
	
	output:
		tuple val(sampleName), path('pile.up') into pile_up_a, pile_up_b
	
	shell:
	"""
		samtools mpileup -aa -A -d 10000 -B -Q 0 --reference $params.referenceSequence -o pile.up resorted.bam
	"""
}



process variantCalling {
	input:
		tuple val(sampleName), path('pile.up') from pile_up_a
	
	output:
		tuple val(sampleName), path('rawVarCalls.tsv') into ivar_out
		
	shell:
	"""
		cat pile.up | ivar variants -p rawVarCalls -g $projectDir/covidRefSequences/covidGenomeAnnotation-NCBI.gff \
		-r $params.referenceSequence -m 10
	"""
}


process kraken2stdDB {
	memory '70 GB'
	cpus 10
	
	input:
		tuple val(sampleName), file('R1.fastq.gz'), file('R2.fastq.gz') from input_fq_b
	
	output:
		tuple val(sampleName), path('k2-std.out') into k2_std_out
		
	shell:
	"""
		numThreads=\$(nproc)
		if $isPairedEnd; then
			kraken2 --paired R1.fastq.gz R2.fastq.gz --threads \$numThreads --report k2-std.out > /dev/null
		else
			kraken2 R1.fastq.gz --threads \$numThreads --report k2-std.out > /dev/null		
		fi	
	"""
}


// The pileup file is parsed to calculate the positionwise quality and depth parameters.
// The result is stored as png files that are added to the html report
process plotCoverageQC {
	input:
		tuple val(sampleName), path('pile.up') from pile_up_b
	
	output:
		tuple val(sampleName), path('pos-coverage-quality.tsv'), path('coverage.png'), path('depthHistogram.png'), path('quality.png'), path('qualityHistogram.png'), path('terminiDensity.png') into QChists
		
	shell:
	"""
		primerBedFileLocal=$primerBedFile
		gapfilename=\${primerBedFileLocal%.*}.uncovered
		$projectDir/plotCoverageQualityPerPos.py pile.up ./ \$gapfilename
	"""
}


// Draw a histogram of all read lengths
process readLengthHist {
	input:
		tuple val(sampleName), file('R1.fastq.gz'), file('R2.fastq.gz') from input_fq_d
	
	output:
		tuple val(sampleName), path('readLengthHist.png') into readLengthHist_png
		
	shell:
	"""
		gzip -dc R1.fastq.gz > allreads.fastq
		if $isPairedEnd; then
			gzip -dc R2.fastq.gz >> allreads.fastq
		fi
		cat allreads.fastq | awk 'NR%4==2' | awk "{ print length}" | $projectDir/plotLengthHist.py
		rm allreads.fastq
	"""
}



// ////////////////////////////////////////////
// VARIANT CALLING
// ////////////////////////////////////////////

process krakenVariantCaller {
	cpus 10
	
	input:
		tuple val(sampleName), env(numReads), path('resorted.fastq.gz') from resorted_fastq_gz_a
	
	output:
		tuple val(sampleName), path('k2-allCovid_bracken.out'), path('k2-majorCovid_bracken.out'), path('k2-allCovid.out'), path('k2-majorCovid.out') into k2_covid_out
	
	shell:
	"""
		numThreads=\$(nproc)
		
		# Check the number of reads. Ignore if there are too few reads
		if [[ \$numReads -gt 10 ]]; then
			kraken2 resorted.fastq.gz --db $projectDir/customDBs/allCovidDB --threads \$numThreads --report k2-allCovid.out > /dev/null
			if [[ \$(cat k2-allCovid.out | wc -l) -eq 1 ]]; then
				# There is a bug in our bracken that fails if no hits.
				echo 100.00\$'\t'0\$'\t'0\$'\t'R\$'\t'1\$'\t'root > k2-allCovid_bracken.out
			else
				bracken -d $projectDir/customDBs/allCovidDB -i k2-allCovid.out -o allCovid.bracken -l P
			fi

			kraken2 resorted.fastq.gz --db $projectDir/customDBs/majorCovidDB --threads \$numThreads --report k2-majorCovid.out > /dev/null
			if [[ \$(cat k2-allCovid.out | wc -l) -eq 1 ]]; then
				# There is a bug in our bracken that fails if no hits.
				echo 100.00\$'\t'0\$'\t'0\$'\t'R\$'\t'1\$'\t'root > k2-majorCovid_bracken.out
			else
				bracken -d $projectDir/customDBs/majorCovidDB -i k2-majorCovid.out -o majorCovid.bracken -l P
			fi
		else
			echo 100.00\$'\t'0\$'\t'0\$'\t'R\$'\t'1\$'\t'root > k2-allCovid_bracken.out
			cp k2-allCovid_bracken.out k2-majorCovid_bracken.out
			cp k2-allCovid_bracken.out k2-allCovid.out
			cp k2-allCovid_bracken.out k2-majorCovid.out
		fi
	"""
}



process pangolinVariantCaller {
	cpus 4
	
	input:
		tuple val(sampleName), env(numReads), path('resorted.bam') from resorted_bam_c
	
	output:
		tuple val(sampleName), env(consensusLineage), path('lineage_report.csv') into pangolin_out
	
	shell:
	"""
		numThreads=\$(nproc)
		
		bcftools mpileup -d 10000 -Ou -f $params.referenceSequence resorted.bam | bcftools call --ploidy 1 -mv -Oz -o calls.vcf.gz
		bcftools index calls.vcf.gz
		cat $params.referenceSequence | bcftools consensus calls.vcf.gz > consensus.fa
		
		pangolin --alignment consensus.fa --threads \$numThreads --outdir ./
		
		# Characterisation of the consensus sequence based on Pangolin output
		# Calculation of the consensus sequence is used to determine the predominant lineage.
		# If available, use the WHO label.
		consensusLineage=\$(tail -n 1 lineage_report.csv | awk -F "," '{ print \$2 }')
		WHOlabel=\$(cat $projectDir/pangolin2WHOlabel.txt | grep \$consensusLineage | awk -F " " '{ print \$2 }')
		if  [[ -n \$WHOlabel ]]; then
			consensusLineage=\$WHOlabel
		fi
	"""
}


process linearDeconVariantCaller {
	cpus 4
	
	input:
		tuple val(sampleName), path('rawVarCalls.tsv') from ivar_out
	
	output:
		tuple val(sampleName), path('linearDeconvolution_abundance.csv'), path('mutationTable.html'), path('VOC-VOIsupportTable.html'), env(mostAbundantVariantPct), env(mostAbundantVariantName), env(linRegressionR2) into linearDeconvolution_out
		
	shell:
	"""
		deconvolutionOutput=\$($projectDir/deconvolveVariants.py rawVarCalls.tsv ./ $params.variantDBfile)
		
		mostAbundantVariantPct=\$(echo \$deconvolutionOutput | awk '{ print \$1 }')
		mostAbundantVariantName=\$(echo \$deconvolutionOutput |awk '{ print \$2 }')
		linRegressionR2=\$(echo \$deconvolutionOutput | awk '{ print \$3 }')
	"""
}


process kallistoVariantCaller {
	cpus 10
	
	input:
		tuple val(sampleName), env(numReads), path('resorted.fastq.gz') from resorted_fastq_gz_b

	output:
		tuple val(sampleName), path('kallisto_abundance.tsv') into kallisto_out
		
	shell:
	"""
		# Check the number of reads. Ignore if there are too few reads
		if [[ $task.attempt -lt 2 ]] && [[ \$numReads -gt 10 ]]; then
			numThreads=\$(nproc)
			kallisto quant --index $projectDir/customDBs/variants.kalIdx --output-dir ./ \
					--plaintext -t \$numThreads --single -l 300 -s 50 resorted.fastq.gz || true
		else
			echo target_id\$'\t'length\$'\t'eff_length\$'\t'est_counts tpm > abundance.tsv
			echo other\$'\t'29903\$'\t'29903\$'\t'0\$'\t'NaN >> abundance.tsv
		fi
		mv abundance.tsv kallisto_abundance.tsv
	"""
}


// https://github.com/rvalieris/LCS
process LCSvariantCaller {	
	input:
		tuple val(sampleName), env(numReads), path('resorted.fastq.gz') from resorted_fastq_gz_c 
	
	output:
		tuple val(sampleName), path('LCS/outputs/decompose/lcs.out') into lcs_out

	shell:
	"""
		if [[ $task.attempt -lt 2 ]] && [[ \$numReads -gt 10 ]]; then
			echo Fetching the LCS repository...
			git clone https://github.com/rvalieris/LCS.git
			cd LCS
			rm .git -rf
			
			echo Preparing the DB...
			mkdir -p outputs/variants_table
			zcat data/pre-generated-marker-tables/pango-designation-markers-v1.2.124.tsv.gz > outputs/variants_table/pango-markers-table.tsv
			
			# A maximum of 100 000 reads are fed into the algorithm to keep the computation time reasonable.
			echo Preparing the sample dataset...
			mkdir data/fastq
			gzip -dc ../resorted.fastq.gz | head -n 400000 > data/fastq/resorted.fastq
			gzip data/fastq/resorted.fastq
			echo "resorted" > data/tags_pool_lcs
			
			echo Executing LCS...
			snakemake --config markers=pango dataset=lcs --cores 2
		else
			echo Not enough covid reads for LCS, skipped.
			mkdir -p LCS/outputs/decompose
			echo sample\$'\t'variant_group\$'\t'proportion\$'\t'mean\$'\t'std_error > LCS/outputs/decompose/lcs.out
			echo ERROR\$'\t'Other\$'\t'1\$'\t'1\$'\t'1 >> LCS/outputs/decompose/lcs.out
		fi
	"""
}


process freyjaVariantCaller {
	input:
		tuple val(sampleName), env(numReads), path('resorted.bam') from resorted_bam_d
		
	output:
		tuple val(sampleName), path('freyja.demix') into freyja_out
		
	// Due to a potential bug, some big fastqs result in a pandas error.
	// Start by generating an empty file to circumvent such failure cases
	shell:
	"""
		if [[ $task.attempt -lt 2 ]] && [[ \$numReads -gt 100 ]]; then
			echo Pileup generation for Freyja
			freyja variants resorted.bam --variants freyja.variants.tsv --depths freyja.depths --ref $params.referenceSequence
			
			echo Demixing variants by Freyja
			freyja demix freyja.variants.tsv freyja.depths --output freyja.demix
		else
			echo FATAL ERROR > freyja.demix
			echo summarized\$'\t'"[('Other', 1.00)]" >> freyja.demix
		fi
	"""
}



// A metadata fetch attemp from NCBI via Entrez-Direct
// Only works if file names explicitly carries an SRR number.
process getNCBImetadata {
	// Allow access to this section only 1 thread at a time to avoid network congestion.
	executor = 'local'
	queueSize = 1
	submitRateLimit = '3/1min'
	
	// NCBI bandwidth limit might cause lookup failures. If so, the next attempt should start with a time delay.
	// Wait some random time so that threads go out of sync.
	errorStrategy { sleep(1000 * Math.random() as long); return 'retry' }
	maxRetries = 55
	
	input:
		tuple val(sampleName), file('R1.fastq.gz'), file('R2.fastq.gz') from input_fq_c
		
	output:
		tuple val(sampleName), env(libraryProtocol), env(seqInstrument), env(isolate), env(collectionDate), env(collectedBy), env(sequencedBy), env(sampleLatitude), env(sampleLongitude), env(sampleLocation) into metadata
	
	shell:
	"""
	srrNumber=$sampleName
	if [[ $task.attempt -lt 50 ]] && [[ \${srrNumber:0:3} == 'SRR' ]]; then
		# The tool returns error: too many requests, bypassing by redirection of error		
		sraQueryResult=\$(esearch -db sra -query \$srrNumber 2>/dev/null)
		sleep 1
		
		if echo \$sraQueryResult | grep -q "<Count>1</Count>"; then
			# Get runinfo from SRA
			echo Downloading metadata for \$srrNumber...	
			echo "\$sraQueryResult" | efetch --format runinfo
			sleep 1
			SRRmetadata=\$(echo "\$sraQueryResult" | efetch --format runinfo 2>/dev/null | grep \$srrNumber)
			
			echo Parsing...
			libraryProtocol=\$(echo \$SRRmetadata | awk -F ',' '{print \$13}')
			seqInstrument=\$(echo \$SRRmetadata | awk -F ',' '{print \$20}')
			isolate=\$(echo \$SRRmetadata  | awk -F ',' '{print \$30}')	
			
			# Get metadata out of biosample db
			echo Fetching biosample data...
			SAMN=\$(echo \$SRRmetadata | awk -F ',' '{print \$26}')
			SAMNmetadata=\$(efetch -db biosample -id \$SAMN 2>/dev/null)
			
			echo Parsing...
			collectionDate=\$(echo "\$SAMNmetadata" | grep "collection date" | awk -F '"' '{print \$2}')
			collectedBy=\$(echo "\$SAMNmetadata" | grep "collected by" | awk -F '"' '{print \$2}')
			sequencedBy=\$(echo "\$SAMNmetadata" | grep SEQUENCED_BY | awk '{ \$1=""; print \$0 }')
			sampleLatitude=\$(echo "\$SAMNmetadata" | grep "latitude and longitude" | awk -F '"' '{print \$2}'\
								| awk '{ print \$1\$2 }')
			sampleLongitude=\$(echo "\$SAMNmetadata" | grep "latitude and longitude" | awk -F '"' '{print \$2}'\
								| awk '{ print \$3\$4 }')
			sampleLocation=\$(echo "\$SAMNmetadata" | grep "geographic location" | awk -F '"' '{print \$2}')
		fi
	fi
	
	if [[ -z \$SAMN ]]; then
		SAMN=Missing
	fi

	if [[ -z \$libraryProtocol ]]; then
		libraryProtocol=Missing
	fi

	if [[ -z \$seqInstrument ]]; then
		seqInstrument=Missing
	fi

	if [[ -z \$isolate ]]; then
		isolate=Missing
	fi

	if [[ -z \$collectionDate ]]; then
		collectionDate=Missing
	fi

	if [[ -z \$collectedBy ]]; then
		collectedBy=Missing
	fi

	if [[ -z \$sequencedBy ]]; then
		sequencedBy=Missing
	fi

	if [[ -z \$sampleLatitude ]]; then
		sampleLatitude="?"
	fi

	if [[ -z \$sampleLongitude ]]; then
		sampleLongitude="?"
	fi

	if [[ -z \$sampleLocation ]]; then
		sampleLocation=Missing
	fi	
	"""
}



// Computation is now mostly over. All threads need to synchronise here.
// We will group based on the sample name and pass everything to the report
// generation steps.
reportInputCh = metadata.join(samtools_stats).join(k2_std_out).join(QChists).join(readLengthHist_png).join(linearDeconvolution_out)
						.join(k2_covid_out).join(pangolin_out).join(kallisto_out).join(freyja_out).join(lcs_out)


///////////////////////////////////////////////
// Report generation and final output
///////////////////////////////////////////////

// Generates a report based on the computation results generated by executeAnalysis.sh
// One separate html report per each sample (i.e. per fastq/fastq pair)
process generateReport {
	input:
		tuple val(sampleName), env(libraryProtocol), env(seqInstrument), env(isolate), env(collectionDate), env(collectedBy), env(sequencedBy), env(sampleLatitude), env(sampleLongitude), env(sampleLocation),
		path('sorted.stats'), path('resorted.stats'),
		path('k2-std.out'),
		path('pos-coverage-quality.tsv'), path('coverage.png'), path('depthHistogram.png'), path('quality.png'), path('qualityHistogram.png'), path('terminiDensity.png'),
		path('readLengthHist.png'),
		path('linearDeconvolution_abundance.csv'), path('mutationTable.html'), path('VOC-VOIsupportTable.html'), env(mostAbundantVariantPct), env(mostAbundantVariantName), env(linRegressionR2),
		path('k2-allCovid_bracken.out'), path('k2-majorCovid_bracken.out'), path('k2-allCovid.out'), path('k2-majorCovid.out'),
		env(consensusLineage), path('lineage_report.csv'),
		path('kallisto_abundance.tsv'),
		path('freyja.demix'),
		path('lcs.out') from reportInputCh
	
	output:
		file "outfolder" into reportCh
	
	shell:
	"""
		export PYTHONHASHSEED=0
		$projectDir/plotPieChartsforAbundance.py ./ $params.variantDBfile linearDeconvolution_abundance.csv \
				kallisto_abundance.tsv k2-allCovid_bracken.out k2-majorCovid_bracken.out freyja.demix lcs.out
		
		export kallistoTopName=\$(cat kallisto.out | sort -k 2 -n | tail -n 1 | awk '{ print \$1 }')
		
		$projectDir/generateReport.sh $sampleName $projectDir/htmlHeader.html $isPairedEnd $params.primerBedFile $projectDir
	"""
}



// The below process runs once per folder and generates a concise summary of all samples after all other
// executions are over.
process summaryPage {
	executor = 'local'
	publishDir "$params.out", mode: 'copy', overwrite: true

	input:
		file 'report' from reportCh.collect()
	
	output:
		file "analysisResults" into analysisResults
	
	shell:
	"""
		$projectDir/generateSummary.sh $projectDir/htmlHeader.html $params.variantDBfile $projectDir
	"""
}


