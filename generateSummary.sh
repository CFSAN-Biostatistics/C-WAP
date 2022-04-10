#!/bin/bash -e

headerFile=$1
variantDBfile=$2
rootDir=$3

summaryfile=./summary.html
cp $headerFile $summaryfile


# This-preprocessing step is necessary to ensure that the results will be reported alphanumerically sorted
# whereas Nextflow reports them in the order of job completion by default.
echo Renaming the output files nicely...
for subdir in $(ls */ -d | tr -d '/'); do
	reportFileName=$subdir/report/report.html
	# Check if this report file exists, else, its contents were already relabelled before
	if [[ -f $reportFileName ]]; then
		# Rename the files/folders such that sample name is a prefix
		sampleName=$(cat $reportFileName | grep "Sample name:" | awk -F'>' '{ print $4 }' | awk -F'<' '{print $1}')
		for file in $(ls $subdir); do
			mv $subdir/$file $subdir/${sampleName}_${file}
		done
	else
		sampleName=$(ls $subdir | head -n 1 | awk -F '_' '{print $1}')
	fi	
	
	mv $subdir $sampleName
done


#######################################################
# Generate a table summarising amplification specificity and the completeness of the genomic coverage.
echo Gathering and tabulating summary statistics...

echo "<h2>Summary</h2>" >> $summaryfile
echo '<table>' >> $summaryfile
echo "<tr>" >> $summaryfile
echo "<th>Sample#</th><th>Sample name</th><th>Total #reads</th><th>Reads aligned PF*</th><th>Genomic coordinates 0X</th><th>Genomic coordinates <10X</th>" >> $summaryfile
echo "</tr>" >> $summaryfile

# Also generate a separate csv file containing the summary table 
csvFile=./summaryTable.csv
echo "Sample#,Sample name,Total #reads,Reads aligned PF,Genomic coordinates 0X,Genomic coordinates <10X, QC_determination, QC_issues" >> $csvFile

sampleNames=()
plottingData=()
coverageFiles=()
for sampleName in $(ls */ -d | tr -d '/'); do
	sampleNames+=($sampleName)
	
	# Extract some of the statistics from individual report files
	reportFileName=$sampleName/${sampleName}_report/report.html
	numReads=$(cat $reportFileName | grep "Total number of reads:" | awk -F '>' '{ print $4 }' | awk -F '<' '{print $1}')
	readsMapped=$(cat $reportFileName | grep "Reads passing filter" | awk -F '>' '{ print $4 }' | awk -F '<' '{ print $1 }')
	pctMapped=$(echo $readsMapped | awk '{print $2}' | tr -d '()%')
	numUncovered=$(cat $reportFileName | grep "All genomic coordinates:" | awk -F '>' '{ print $4 }' | awk -F '<' '{ print $1 }')
	numPoorlyCovered=$(cat $reportFileName | grep "All genomic coordinates:" | awk -F '>' '{ print $6 }' | awk -F '<' '{ print $1 }')
	
	plottingData+=($sampleName $(echo $readsMapped | awk '{print $1}' | tr -d 'nt'))
	coverageFiles+=($sampleName/${sampleName}_pos-coverage-quality.tsv)

	# Record as a row in the data table	
	echo "<tr>" >> $summaryfile
	echo "<td>${#sampleNames[@]}</td><td><a href=\"./$sampleName/${sampleName}_report/report.html\">$sampleName</a></td>" >> $summaryfile 
	echo "<td>$numReads</td><td>$readsMapped</td><td>$numUncovered</td><td>$numPoorlyCovered</td>" >> $summaryfile
	echo "</tr>" >> $summaryfile
	
	QC_determination="not_performed" # not_performed / passed / failed / contains_issues
	QC_flags="None" # None / low_quality_sequence / sequence_contaminated / low_average_coverage / uneven_sequence_coverage
	# low_coverage_breadth / truncated_reads / sequence_amplification_artifacts / signal_to_noise_ratio_low
	
	echo "${#sampleNames[@]},$sampleName,$numReads,$readsMapped,$numUncovered,$numPoorlyCovered,${QC_determination},${QC_flags[@]}" >> $csvFile
done

echo '</table>' >> $summaryfile
echo "<br><br><br>" >> $summaryfile



#######################################################
# Insert a bar plot comparing covid reads in all samples of this run.
echo Generating sample load comparison plots...
$rootDir/plotSampleYieldComparison.py ${plottingData[@]}

echo "<div>" >> $summaryfile
echo "   <div id=\"figdiv\">" >> $summaryfile
echo "      <img src=\"./covidReadsSummary.png\" alt=\"Num. SC2 reads\" width=\"100%\" class=\"center\">" >> $summaryfile
echo "   </div>" >> $summaryfile
echo "   <div id=\"figdiv\" style=\"padding-left: 5%; width: 45%;\">" >> $summaryfile
echo "       *Quantity of raw reads that align to the reference sequence and pass filter, i.e. the read length after adaptor trimming &#8805;30 and minimum read quality &#8805;20 within a sliding window of width 4. SNR refers to the ratio of SC2-mapping reads aligned that pass filter in the sample vs. that in the auto-detected negative control samples (if any). The dashed line represents the baseline level of covid reads detected from the negative control or their average if multiple negative controls we included." >> $summaryfile
echo "   </div>" >> $summaryfile
echo "</div>" >> $summaryfile
echo "<br><br><br>" >> $summaryfile


# Run the coverage stats through the machine learning-based QC model to predict the accuracy of the samples given their coverage pattern.
echo Executing QC-bot...
$rootDir/AI/predictAccuracy.py $rootDir/AI/trainedModel.pkl ${coverageFiles[@]}

echo >> $summaryfile
echo "<h2>QC-bot (Not yet functional)</h2>" >> $summaryfile

echo "<div>" >> $summaryfile
echo "   <div id=\"figdiv\">" >> $summaryfile
echo "      <img src=\"./QC_predictions.png\" alt=\"Num. SC2 reads\" width=\"100%\" class=\"center\">" >> $summaryfile
echo "   </div>" >> $summaryfile
echo "   <div id=\"figdiv\" style=\"padding-left: 5%; width: 45%;\">" >> $summaryfile
echo "       Machine-learning based prediction of the SC2 variant calling accuracy of Freyja of this dataset. The model is a random forest trained on FDA/CFSAN's experimental wastewater WGS data obtained in January 2022 and aims to assess the impact of the potential coverage gaps on the variant abundance estimates. The plotted values represent the predicted deviation of the omicron percentage points from the value that would have been obtained if the coverage was near-complete." >> $summaryfile
echo "   </div>" >> $summaryfile
echo "</div>" >> $summaryfile
echo "<br><br><br>" >> $summaryfile



#######################################################
# Place the quality plots, one figure per sample on a nx2 grid.
echo Adding figures...
for ((i = 0; i < ${#sampleNames[@]}; i+=2)); do
	echo >> $summaryfile

	echo "<span style=\"float:left;\"><a href=\"./${sampleNames[$i]}/${sampleNames[$i]}_report/report.html\">${sampleNames[$i]}</a></span>" >> $summaryfile
	echo "<span style=\"float:right;\"><a href=\"./${sampleNames[$(expr $i + 1)]}/${sampleNames[$(expr $i + 1)]}_report/report.html\">" >> $summaryfile
	echo "${sampleNames[$(expr $i + 1)]}</a></span><br>" >> $summaryfile
	
	echo "<img src=\"./${sampleNames[$i]}/${sampleNames[$i]}_report/quality.png\" width=\"49%\" class=\"center\">" >> $summaryfile
	echo "<img src=\"./${sampleNames[$(expr $i + 1)]}/${sampleNames[$(expr $i + 1)]}_report/quality.png\" width=\"49%\" class=\"center\">" >> $summaryfile
	echo "<br><br><br>" >> $summaryfile
done

	
#######################################################
#echo Appending a list of data analysis parameters...
#echo >> $summaryfile
#echo "<h2>Software configuration</h2>" >> $summaryfile

#pangolinVersion=$(pangolin -v)
#pangolearnVersion=$(pangolin -pv)

#if [[ $platform == "Illumina" ]]; then
#	echo "Bowtie2 v)bowtie2 --version | head -n 1 | awk '{ print $3 }'), " >> $summaryfile
#else
#	echo "Minimap2:" $(minimap2 --version) >> $summaryfile
#fi

#echo $(samtools --version | head -n 2)", " >> $summaryfile
#echo $(ivar version | head -n 1)", " >> $summaryfile
#echo $(kraken2 -v | head -n 1)", " >> $summaryfile
#echo $(kallisto version)"." >> $summaryfile

#allIncludedLineages=$($rootDir/listVariantsAvail.py $variantDBfile)
#echo Lineage definitions were compiled on $(date +%Y-%m-%d -r $variantDBfile) \
#		"from <a href=\"https://github.com/cov-lineages/constellations/tree/main/constellations/definitions\">constellations</a>." >> $summaryfile
#echo Lineage signature file was compiled on $(date +%Y-%m-%d -r $variantDBfile) \
#		and includes lineages: $allIncludedLineages. >> $summaryfile
#echo Lineage assignment to the consensus sequence was performed by $pangolinVersion using the classification tree of $pangolearnVersion. >> $summaryfile


#cp $headerFile $summaryfile
echo "</body>" >> $summaryfile
echo "</html>" >> $summaryfile


#######################################################
echo Consolidating reports...
wkhtmltopdf --enable-local-file-access --page-size Letter --margin-top 10mm --margin-bottom 0 \
		--margin-left 0 --margin-right 0 --print-media-type --title "Wastewater report" \
		summary.html summary.pdf

gs -dNOPAUSE -dPDFSETTINGS=/prepress -sDEVICE=pdfwrite -sOUTPUTFILE=consolidated.pdf \
	-dBATCH summary.pdf ./*/*report/report.pdf


mkdir analysisResults
shopt -s extglob
mv ./!(analysisResults) ./analysisResults/

