#!/bin/bash -e

headerFile=$1
variantDBfile=$2
rootDir=$3

summaryfile=./summary.html
cp $headerFile $summaryfile



# This-preprocessing step is necessary to ensure that the results will be reported alphanumerically sorted
# whereas Nextflow reports them in the order of job completion by default.
echo Renaming the output files nicely...
for subdir in `ls */ -d | tr -d '/'`; do
	reportFileName=$subdir/report/report.html
	sampleName=`cat $reportFileName | grep "Sample name:" | awk -F'>' '{ print $4 }' | awk -F'<' '{print $1}'`

	# Rename the files/folders such that sample name is a prefix
	for file in `ls $subdir`; do
		mv $subdir/$file $subdir/${sampleName}_${file}
	done
	mv $subdir $sampleName
done


#######################################################
# Generate a table summarising amplification specificity and the completeness of the genomic coverage.
echo Gathering and tabulating summary statistics...

echo "<h2>SUMMARY</h2>" >> $summaryfile
echo '<table>' >> $summaryfile
echo "<tr>" >> $summaryfile
echo "<th>Sample#</th><th>Sample name</th><th>Total #reads</th><th>% covid hits</th><th>Genomic coordinates 0X</th><th>Genomic coordinates <10X</th>" >> $summaryfile
echo "</tr>" >> $summaryfile

counter=0
sampleNames=()
for sampleName in `ls */ -d | tr -d '/'`; do
	sampleNames[$counter]=$sampleName
	counter=$((counter+1))
	
	# Extract some of the statistics from individual report files
	reportFileName=$sampleName/${sampleName}_report/report.html
	numReads=`cat $reportFileName | grep "Total number of reads:" | awk -F'>' '{ print $4 }' | awk -F'<' '{print $1}'`
	pctCovid=`cat $reportFileName | grep "Hits to SARS-Cov2 genome" | awk -F'(' '{ print $3 }' | awk -F'%' '{ print $1 }'`
	numUncovered=`cat $reportFileName | grep "ncovered genomic coordinates (0X):" | awk -F'>' '{ print $4 }' | awk '{ print $1 }'`
	numPoorlyCovered=`cat $reportFileName | grep "poorly covered genomic coordinates (<10X):" | awk -F'>' '{ print $4 }' | awk '{ print $1 }'`

	# Record as a row in the data table
	echo "<tr>" >> $summaryfile
	echo "<td>$counter</td><td><a href=\"./$sampleName/${sampleName}_report/report.html\">$sampleName</a></td>" >> $summaryfile 
	echo "<td>$numReads</td><td>$pctCovid</td><td>$numUncovered</td><td>$numPoorlyCovered</td>" >> $summaryfile
	echo "</tr>" >> $summaryfile
done

echo '</table>' >> $summaryfile
echo "<br><br><br>" >> $summaryfile


# Place the quality plots, one figure per sample on a nx2 grid.
echo Adding figures
for ((i = 0; i < ${#sampleNames[@]}; i+=2)); do
	echo "<span style=\"float:left;\"><a href=\"./${sampleNames[$i]}/${sampleNames[$i]}_report/report.html\">${sampleNames[$i]}</a></span>" >> $summaryfile
	echo "<span style=\"float:right;\"><a href=\"./${sampleNames[`expr $i + 1`]}/${sampleNames[`expr $i + 1`]}_report/report.html\">" >> $summaryfile
	echo "${sampleNames[`expr $i + 1`]}</a></span><br>" >> $summaryfile
	
	echo "<img src=\"./${sampleNames[$i]}/${sampleNames[$i]}_report/quality.png\" width=\"49%\" class=\"center\">" >> $summaryfile
	echo "<img src=\"./${sampleNames[`expr $i + 1`]}/${sampleNames[`expr $i + 1`]}_report/quality.png\" width=\"49%\" class=\"center\">" >> $summaryfile
	echo "<br><br><br>" >> $summaryfile
done
	
	
#######################################################
echo Appending a list of data analysis parameters...
echo >> $summaryfile
echo "<h2>Software configuration</h2>" >> $summaryfile

pangolinVersion=`pangolin -v`
pangolearnVersion=`pangolin -pv`

if [[ $platform == "Illumina" ]]; then
	echo "Bowtie2 v`bowtie2 --version | head -n 1 | awk '{ print $3 }'`, " >> $summaryfile
else
	echo "Minimap2:" `minimap2 --version` >> $summaryfile
fi

echo `samtools --version | head -n 2`", " >> $summaryfile
echo `ivar version | head -n 1`", " >> $summaryfile
echo `kraken2 -v | head -n 1`", " >> $summaryfile
echo `kallisto version`"." >> $summaryfile

allIncludedLineages=`$rootDir/listVariantsAvail.py $variantDBfile`
echo "Lineage definitions were compiled on" `date +%Y-%m-%d -r $variantDBfile` \
		"from <a href=\"https://github.com/cov-lineages/constellations/tree/main/constellations/definitions\">constellations</a>." >> $summaryfile
echo "Lineage signature file was compiled on" `date +%Y-%m-%d -r $variantDBfile` \
		"and includes lineages: $allIncludedLineages." >> $summaryfile
echo "Lineage assignment to the consensus sequence was performed by $pangolinVersion using the classification tree of $pangolearnVersion." >> $summaryfile
	
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

