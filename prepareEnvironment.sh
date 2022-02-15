#!/bin/bash

# FDA/CFSAN-specific script
# Prepares the bash environment by loading necessary modules or conda packages
# Please adjust as needed by your system to make sure all executables are in the path
# before invoking the pipeline.


# Exit if any of these commands fail
set -e


echo Loading necessary modules...
module purge
module load /nfs/software/modules/kraken2/2.1.2 bracken/2.5
module load samtools/1.13 ivar
module load bowtie2
module load minimap2/2.22
module load edirect
module load nextflow


echo Loading conda modules...
condaBin=/projects/covidtrakr/software/miniconda3/bin/
export PATH=$PATH:$condaBin

pangolinDir=/projects/covidtrakr/software/miniconda3/envs/pangolin/bin
export PATH=$PATH:$pangolinDir


# Freyja's python library is incompatible with the libraries on raven2
# so everything to be used out of its own conda package.
freyja () {
	originalPath=$PYTHONPATH
	unset PYTHONPATH
	/projects/covidtrakr/software/miniconda3/envs/freyja-env/bin/freyja $@
	export PYTHONPATH=$originalPath
}
export -f freyja


if ! command -v wkhtmltopdf; then
	echo ERROR: wkhtmltopdf does not exist.
	echo Hint: You can download from: "https://wkhtmltopdf.org/downloads.html"
	exit 1
fi


export organisation=fda.hhs.gov


# Unset "stop on exit" behaviour so that this script can be sourced for other purposes.
set +e

echo The environment is set

