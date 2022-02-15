#! /bin/bash	

# Prepare an html report. The html file header and standard page banner are copied from a stored file for simplicity.


sampleName=$1
headerFile=$2
isPairedEnd=$3

reportFile=report.html
cp $headerFile $reportFile


IFS='.' read -ra username <<< "$USER"
fullName="${username[0]} ${username[1]}"
timestamp=`date +%Y-%m-%d,\\ %T\\ %Z`

echo "<table>" >> $reportFile
echo "<tr>" >> $reportFile
echo "    <td>Sample name:</td><td>$sampleName</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Date generated:</td><td>$timestamp</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Executed by:</td><td>$fullName (<a href=\\"mailto:$USER@fda.hhs.gov?subject=Wastewater report generated on $timestamp\\">$USER@fda.hhs.gov</a>)</td>" >> $reportFile
echo "</tr>" >> $reportFile

#echo "<tr>" >> $reportFile
#echo "    <td>Executed on:</td><td>`hostname -I | awk '{print $1}'` (aka `hostname`)</td>" >> $reportFile
#echo "</tr>" >> $reportFile
echo "</table>" >> $reportFile


#######################################################
echo Compiling read statistics
numAligned=`cat sorted.stats | grep "reads mapped:" | awk '{ print $NF }'`
avgQuality=`cat sorted.stats | grep "average quality" | awk '{ print $NF }'`
avgLength=`cat sorted.stats | grep "average length" | awk '{ print $NF }'`

numPassedQuality=`cat resorted.stats | grep "raw total sequences" | awk '{ print $4 }'`
avgQualityPassed=`cat resorted.stats | grep "average quality" | awk '{ print $NF }'`
avgLengthPassed=`cat resorted.stats | grep "average length" | awk '{ print $NF }'`
let "avgCoveragePassed = $avgLengthPassed * $numPassedQuality / 29903" || true


# Deduce the total number of reads from the kraken2 output
# Use taxid column, i.e. column 5
numUnclassified=`head k2-std.out | awk '$5==0' | awk '{ print $2 }'`
if [[ -z $numUnclassified ]]; then
	numUnclassified=0
fi

numClassified=`head k2-std.out | awk '$5==1' | awk '{ print $2 }'`
if [[ -z $numClassified ]]; then
	numClassified=0
fi

# Calculation of total number of reads
let "numReads = $numUnclassified + $numClassified"
if $isPairedEnd; then
# Correction for undercounting by kraken2 in the paired mode
	let "numReads = 2 * $numReads"
fi


####################################################################
echo >> $reportFile
echo "<br>" >> $reportFile
echo "<h2>Sequencing summary</h2>" >> $reportFile
echo '<table>' >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Sequencing chemistry:</td><td>$libraryProtocol with $seqInstrument</td>" >> $reportFile
echo "</tr>" >> $reportFile

# Catchment site the sample was collected from
echo "<tr>" >> $reportFile
echo "    <td>Source site:</td><td><a href=\"https://www.openstreetmap.org/#map=3/$sampleLatitude/$sampleLongitude\">$sampleLocation\
		($sampleLatitude,$sampleLongitude)</a></td>" >> $reportFile
echo "</tr>" >> $reportFile	

# Date the WW sample was collected from the sanitary system
echo "<tr>" >> $reportFile
echo "    <td>Sampling date:</td><td>$collectionDate</a></td>" >> $reportFile
echo "</tr>" >> $reportFile

# Agency collecting the WW sample
echo "<tr>" >> $reportFile
echo "    <td>Collected by:</td><td>$collectedBy</td>" >> $reportFile
echo "</tr>" >> $reportFile

# Lab that sequenced the sample
echo "<tr>" >> $reportFile
echo "    <td>Sequenced by:</td><td>$sequencedBy</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Total number of reads:</td><td>$numReads</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Reads aligned:</td><td>$numAligned (`expr 100 \\* $numAligned / $numReads`%)</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Average read quality:</td><td>$avgQuality</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Average read length:</td><td>$avgLength</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Reads passing filter:</td><td>$numPassedQuality (`expr 100 \\* $numPassedQuality / $numReads`%)</td>" \
		>> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Average read quality passing filter:</td><td>$avgQualityPassed</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Average read length passing filter:</td><td>$avgLengthPassed</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Average coverage passing filter:</td><td>${avgCoveragePassed}X</td>" >> $reportFile
echo "</tr>" >> $reportFile
echo "</table>" >> $reportFile

echo "<br>" >> $reportFile
echo "A read passes filter if the read length after adaptor trimming &#8805; 30 and minimum \
		read quality &#8805; 20 within a sliding window of width 4." >> $reportFile


#######################################################
echo >> $reportFile
echo "<br><br>" >> $reportFile
echo "<h2>Overall sequence characteristics</h2>" >> $reportFile

echo "<img src=\"./coverage.png\" alt=\"Coverage vs. genome coordinate plot\" width=\"49%\" class=\"center\">" >> $reportFile
echo "<img src=\"./quality.png\" alt=\"Quality vs. genome coordinate plot\" width=\"49%\" class=\"center\">" >> $reportFile
echo "<p> NOTE: The red shaded areas marked with a (<font color="red">*</font>) are not covered by the design of the library preparation kit and hence excluded from analyses. <font color="magenta">Magenta</font> curves represent moving average with a window width of 1kb." >> $reportFile
echo "<br><br><br><br>" >> $reportFile
echo "<img src=\"./depthHistogram.png\" alt=\"Coverage histogram\" width=\"49%\" class=\"center\">" >> $reportFile
echo "<img src=\"./qualityHistogram.png\" alt=\"Quality histogram\" width=\"49%\" class=\"center\">" >> $reportFile
echo "<br><br><br><br>" >> $reportFile
echo "<img src=\"./readLengthHist.png\" alt=\"Read length histogram\" width=\"49%\" class=\"center\">" >> $reportFile
echo "<img src=\"./terminiDensity.png\" alt=\"Read termini density\" width=\"49%\" class=\"center\">" >> $reportFile

echo >> $reportFile
echo "<br>" >> $reportFile

# Check for coverage data. If the coverage is too low, any variant calling will be inaccurate
if [ $avgCoveragePassed -lt 10 ]; then
	echo "<p> <font color=\"red\"> WARNING: The sequence coverage is very low (${avgCoveragePassed}X) </font>" >> $reportFile
fi

numUncoveredLociByDesign=0
primerBedFileLocal=$params.primerBedFile
gapfilename=${primerBedFileLocal%.*}.uncovered
while IFS= read -r line; do
	gapBegin=`echo $line | awk '{ print $1 }'`
	gapEnd=`echo $line | awk '{ print $2 }'`
	let "numUncoveredLociByDesign+=$gapEnd - $gapBegin + 1"
done < $gapfilename
numUncoveredLoci=`cat pos-coverage-quality.tsv | awk '{ print $2 }' | grep -w 0 | wc -l`
numPoorlyCoveredLoci=`cat pos-coverage-quality.tsv | awk '{ print $2 }' |
							grep -w -e 0 -e 1 -e 2 -e 3 -e 4 -e 5 -e 6 -e 7 -e 8 -e 9 | wc -l`

numCovid=`cat k2-std.out | grep Orthocoronavirinae | grep -v unclass | awk '{ print $2 }'`
numHuman=`cat k2-std.out | grep sapiens | awk '{ print $2 }'`
numSynthetic=`cat k2-std.out | grep other | grep 28384 | awk '{ print $2 }'`

if [[ -z $numCovid ]]; then
	pctCovid=0.00
	numCovid=0
else
	pctCovid=`cat k2-std.out | grep Orthocoronavirinae | grep -v unclass | awk '{ print $1 }'`
fi

if [[ -z $numHuman ]]; then
	pctHuman=0.00
	numHuman=0
else
	pctHuman=`cat k2-std.out | grep sapiens | awk '{ print $1 }'`
fi

if [[ -z $numSynthetic ]]; then
	pctSynthetic=0.00
	numSynthetic=0
else
	pctSynthetic=`cat k2-std.out | grep other | grep 28384 | awk '{ print $1 }'`
fi


echo '<table>' >> $reportFile
echo "<tr>" >> $reportFile
echo "    <td># uncovered genomic coordinates (0X) by kit design:</td><td>$numUncoveredLociByDesign nt 
		(`expr 100 \\* $numUncoveredLociByDesign / 29903`%)</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td># uncovered genomic coordinates (0X):</td><td>$numUncoveredLoci nt (`expr 100 \\* $numUncoveredLoci / 29903`%)</td>" 			>> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td># poorly covered genomic coordinates (<10X):</td><td>$numPoorlyCoveredLoci nt (`expr 100 \\* $numPoorlyCoveredLoci / 29903`%)</td>" \
		>> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Hits to SARS-Cov2 genome (kraken2):</td><td>$numCovid reads (${pctCovid}%) </td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Hits to human genome (kraken2):</td><td>$numHuman reads (${pctHuman}%)</td>" >> $reportFile
echo "</tr>" >> $reportFile

echo "<tr>" >> $reportFile
echo "    <td>Hits to synthetic sequences (kraken2, taxid 28384):</td><td>$numSynthetic reads (${pctSynthetic}%)</td>" >> $reportFile
echo "</tr>" >> $reportFile
echo "</table>" >> $reportFile


echo "<br>" >> $reportFile
echo "<h2>Detected variants (Experimental)</h2>" >> $reportFile
echo "<div>" >> $reportFile
echo "    <div id=\"figdiv\">" >> $reportFile
echo "        <img src=\"./pieChart_deconvolution.png\" alt=\"Abundance of variants by deconvolution\" width=\"100%\">" >> $reportFile
echo "    </div>" >> $reportFile
echo "    <div id=\"figdiv\">" >> $reportFile
echo "        <img src=\"./pieChart_kallisto.png\" alt=\"Abundance of variants by kallisto\" width=\"100%\">" >> $reportFile
echo "    </div>" >> $reportFile
echo "</div>" >> $reportFile

echo "<div id=\"caption\">" >> $reportFile
echo "    <p> Based on deconvolution, <a href=\"https://outbreak.info/situation-reports?pango=$mostAbundantVariantName\">$mostAbundantVariantName</a> 
	is estimated to constitute $mostAbundantVariantPct% of the viral\
	particles and hence is the most abundant variant in the sample. The R<sup>2</sup> for the linear regression was $linRegressionR2.\
	Variants that were detected less than 5% were grouped under \"Other\"" >> $reportFile
echo "    <p> Based on the consensus sequence of the observed reads, the \"ensemble-averaged sequence\" most closely resembles \
	the <a href=\"https://outbreak.info/situation-reports?pango=$consensusLineage\">$consensusLineage</a> lineage. \
	If this is a sample consisting of a single source of pathogens or an overwhelming majority of the different sources \
	are infected with the same variant, the sample is dominated by this variant." >> $reportFile
echo "    <p> Based on mapping individual reads to the variant consensus sequences in the reference database, \
	kallisto predicts that the sample is dominated by <a href=\"https://outbreak.info/situation-reports?pango=$kallistoTopName\"\
	>$kallistoTopName</a> lineage. Accuracy of this measure is expected to improve if the \
	input data consists of long reads as opposed to convolution." >> $reportFile
echo "</div>" >> $reportFile
echo "<br>" >> $reportFile
	
echo "<br>" >> $reportFile
echo "<div>" >> $reportFile
echo "    <div id=\"figdiv\">" >> $reportFile
echo "        <img src=\"./pieChart_k2_allCovid.png\" alt=\"Abundance of variants by kraken2+bracken\"\
					width=\"100%\">" >> $reportFile
echo "    </div>" >> $reportFile
echo "    <div id=\"figdiv\">" >> $reportFile
echo "        <img src=\"./pieChart_k2_majorCovid.png\" alt=\"Abundance of variants by kraken2+bracken\"\
					width=\"100%\">" >> $reportFile
echo "    </div>" >> $reportFile
echo "</div>" >> $reportFile
echo "<br>" >> $reportFile

echo "<br>" >> $reportFile
echo "<div>" >> $reportFile
echo "    <div id=\"figdiv\">" >> $reportFile
echo "        <img src=\"./pieChart_freyja.png\" alt=\"Abundance of variants by Freyja\"\
				width=\"100%\">" >> $reportFile
echo "    </div>" >> $reportFile
echo "</div>" >> $reportFile
echo "<br>" >> $reportFile

# Append the VOC - VOI support table generated by the above Python script to the report
cat VOC-VOIsupportTable.html >> $reportFile
rm VOC-VOIsupportTable.html


#######################################################
echo Appending a detailed list of all detected mutations...
echo >> $reportFile
echo "<h2>Detected mutations</h2>" >> $reportFile
echo "Only genomic coordinates with at least 10X coverage were considered." >> $reportFile
echo "<br>" >> $reportFile

echo "<table>" >> $reportFile
echo "<tr>" >> $reportFile
echo "    <th>Position</th> <th>Ref. base</th> <th>Alt. base</th> <th>Alt. freq</th> <th>p-value</th> <th>Mutation name</th> <th>Compatible lineages</th>" >> $reportFile
echo "</tr>" >> $reportFile
	
# Copy the list of mutations from the temporary file generated by python
cat mutationTable.html >> $reportFile
rm mutationTable.html
echo "</table>" >> $reportFile
echo "<br>" >> $reportFile
	
echo >> $reportFile
echo "</body>" >> $reportFile
echo "</html>" >> $reportFile


######################################
echo Generating a pdf version of the report...
wkhtmltopdf --enable-local-file-access --page-size Letter --margin-top 10mm --margin-bottom 0 \
	--margin-left 0 --margin-right 0 --print-media-type --title "Wastewater report" \
	report.html report.pdf


# Re-organise the folder so that it is ready to export.
mkdir report
mv report.* report/
mv *.png report/

mkdir outfolder
shopt -s extglob
mv ./!(outfolder) ./outfolder/

