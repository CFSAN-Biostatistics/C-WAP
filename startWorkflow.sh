#!/bin/bash

# Covid Wastewater Analysis Code
# A driver script that parses the input arguments, and starts the analysis and reporting codes

# Exit if any of these commands fail
set -e


usageshort() {
	echo "usage: startWorkflow.sh [-h] [-qgpn] [-f FILE] [-r FILE] [-s FILE] [-b FILE] [-w FILE] [-o DIR]"
	echo '  -h for detailed help message'
	echo
}



usagelong() {
	echo "usage: startWorkflow.sh [-h] [-qpn] [-f FILE] [-r FILE] [-s FILE] [-b FILE] [-w FILE] [-o DIR]"
	echo
	echo "Example usage: "
	echo "    startWorkflow.sh -f sample1_R1.txt -r sample2_R2.fastq -o ./outputDir"
	echo
	echo 'Run the covid19 workflow that will import the fastq files, trim adaptors, apply quality trimming, call variants and 
			generate a pdf report.'
	echo
	echo
	echo 'Options:'
	echo '  -h             : Show this help message and exit.'
	echo
	echo '  -q             : Perform QC only, do not estimate variant abundance.'
	echo
	echo '  -g             : Perform the analysis task as a job on the computation grid.'
	echo
	echo '  -p	           : Use parameters optimised for PacBio platform'
	echo
	echo '  -n	           : Use parameters optimised for the ONT platform'
	echo
	echo '  -f FILE        : Full path to the fastq file containing the forward sequencing reads obtained on an Illumina platform (R1).'
	echo
	echo '  -r FILE        : Full path to the fastq file containing the reverse sequencing reads obtained on an Illumina platform (R2).'
	echo
	echo '  -s FILE        : Full path to the single fastq file containing either interleaved Illumina reads or'
	echo '                         long reads obtained on Pacbio or Oxford Nanopore platforms.'
	echo
	echo '  -b FILE        : Full path to the bed file listing the primers.'	
	echo
	echo "  -w FILE        : Full path to the reference sequence\'s fasta file (default=Wuhan)."
	echo
	echo '  -o STR         : Output directory for all interim results and all other output files.'
	echo
}


unset R1filename R2filename singlefilename primerBedFile


# --------------------------------------------------------
# getopts command line option handler:
# ---------------------------------------------------------
platform=Illumina
performQConly=false
useHPC=false

while getopts ":hqgpnf:r:s:o:w:b:" option; do
	if [ "$option" == ":" ]; then
	echo "Missing argument for option -- '$OPTARG'"
		echo
		usageshort
		exit 1
	fi

	case $option in
		q)
			performQConly=true
			;;
		g)
			useHPC=true
			;;
		p)
			if ! [[ $platform == "Illumina" ]]; then
				echo Conflicting platform information provided.
				exit 1
			else
				platform=Pacbio
			fi
			;;
		n)
			if ! [[ $platform == "Illumina" ]]; then
				echo Conflicting platform information provided.
				exit 1
			else
				platform=ONT
			fi
			;;
		h)
			usagelong
			exit 0
			;;
		\?)
			echo "Invalid option -- '$OPTARG'"
			echo
			usageshort
			exit 1
			;;
		":")
			echo "Missing argument for option -- '$OPTARG'"
			echo
			usageshort
			exit 2
			;;
		f)
			if ! [[ -f  "$OPTARG" ]]; then
				echo File $OPTARG does not exist.
				exit 1
			else
				R1filename=`readlink -f $OPTARG`
			fi
			;;
		r)
			if ! [[ -f  "$OPTARG" ]]; then
				echo File $OPTARG does not exist.
				exit 1
			else
				R2filename=`readlink -f $OPTARG`
			fi
			;;
		s)
			if ! [[ -f  "$OPTARG" ]]; then
				echo File $OPTARG does not exist.
				exit 1
			else
				singlefilename=`readlink -f $OPTARG`
			fi
			;;
		b)
			if ! [[ -f  "$OPTARG" ]]; then
				echo "File $OPTARG does not exist."
				exit 1
			else
				primerBedFile=`readlink -f $OPTARG`
			fi
			;;
		o)
			outDir=`readlink -f $OPTARG`
			if ! [[ -d "$outDir" ]]; then
				echo "Output folder $outDir does not exist. Created."
				mkdir $outDir
			else
				echo "Output folder $outDir already exists, the contents will potentially be overwritten."
			#	rm -r $outDir
			#	mkdir $outDir
			fi
			;;
		w)
			if ! [[ -f  "$OPTARG" ]]; then
				echo "File $OPTARG does not exist."
				exit 1
			else
				referenceSequence=`readlink -f $OPTARG`
			fi
			;;
	esac
done



if [ -z "$outDir" ]; then
    echo "No output directory has been provided. Defaulting to ./out"
	outDir=`readlink -f ./out`
fi


# Path to the folder containing the scripts
path2scripts=`dirname $0`
cd $path2scripts


# ---------------------------------------------------------
# Validation for command line arguments
# ---------------------------------------------------------

# A valid single fastq or Illumina read pair needs to exist
if  [[ -n $singlefilename ]] && [[ -n $R1filename || -n $R2filename ]]; then
	echo
	echo "Too many input parameters (single/interleaved or paired?)."
	echo
	usageshort
	exit 1
fi


if [[ $platform != Illumina && -z $singlefilename ]] || [[ $platform == Illumina && -z $singlefilename && -z $R1filename ]]; then
   echo
   echo No/insufficient input to process was provided.
   echo
   usageshort
   exit 1
fi


if  [[ -n $R1filename && -z $R2filename ]]; then
	echo Only R1 file was provided for a paired-end setup. Assumed same file basename and extension.
	R2filename=`echo $R1filename | sed 's/R1/R2/'`
	if ! [[ -f $R2filename ]]; then
		echo
		echo $R2filename does not exist.
		exit 1
	else
		echo "    R2=$R2filename"
	fi
fi



if [[ -z $primerBedFile ]]; then
	case $platform in
		Illumina)
			echo A primer.bed file was not provided, defaulting to NEB VarSkip Short kit.
			primerBedFile=./covidRefSequences/varskipShort.bed
			;;
		Pacbio)
			echo A primer.bed file was not provided, defaulting to NEB VarSkip Long protocol.
			primerBedFile=./covidRefSequences/varskipLong.bed
			;;
		ONT)
			echo A primer.bed file was not provided, defaulting to Freed protocol.
			primerBedFile=./covidRefSequences/freed.bed
			;;
	esac
fi


if [ -z "$referenceSequence" ]; then
    echo "No reference sequence has been provided. Defaulting to Wuhan (NC_045512.2)"
	referenceSequence=./covidRefSequences/wuhan.fa
fi


# Output sample name to the report file
if [[ -n "$singlefilename" ]]; then
	sampleName=`basename $singlefilename | awk -F '_' '{ print $1 }' | awk -F '.' '{ print $1 }'`
else
	sampleName=`basename $R1filename | awk -F '_' '{ print $1 }' | awk -F '.' '{ print $1 }'`
fi



#######################################################
##### Execute the analysis workflow####################

# Set the environment by loading necessary modules, conda packages etc.
# This is system specific, please adjust the provided script according to your own setup.
source ./prepareEnvironment.sh


# Check whether the system is an HPC controlled via SLURM
if $useHPC; then
	# Pangolin and some other 3rd party software use space in /tmp by default,
	# which causes an out of space error if someone else sharing the node uses heavily.
	# If such an issue occurs, uncomment below:
	# export TMPDIR=$outDir

	echo Executing in SLURM mode, submitting job
	jobID=$(sbatch --parsable -N 1 -c 20 --mem-per-cpu=6G --time 2:00:00 -o $outDir/execution.log \
		./executeAnalysis.sh $outDir $referenceSequence $primerBedFile $performQConly \
		$platform $sampleName $R1filename $R2filename $singlefilename)
	
	jobID2=$(sbatch --parsable -N 1 -c 2 --time 0:10:00 --dependency=afterok:$jobID -o $outDir/reporting.log \
		./generateReport.sh $outDir $primerBedFile $performQConly $platform $sampleName \
		$R1filename $R2filename $singlefilename)

	sbatch -N 1 -c 2 --time 0:10:00 --dependency=afterok:$jobID2 -o $outDir/printing.log \
		./html2pdf.sh $outDir $sampleName
else
	# Ordinary desktop (or another scheduler). Execute the job as is.
	echo Executing in desktop mode
	./executeAnalysis.sh $outDir $referenceSequence $primerBedFile $performQConly $platform $sampleName \
		$R1filename $R2filename $singlefilename | tee $outDir/execution.log
	./generateReport.sh $outDir $primerBedFile $performQConly $platform $sampleName $R1filename $R2filename \
		$singlefilename | tee $outDir/reporting.log
	./html2pdf.sh $outDir $sampleName
fi


echo "(Well) done!"
exit 0

