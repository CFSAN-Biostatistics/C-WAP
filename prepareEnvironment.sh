#!/bin/bash

# CFSAN-specific
# Prepares the bash environment by loading necessary modules or conda packages


echo Loading necessary modules...
module load /nfs/software/modules/kraken2/2.1.2 bracken
module load samtools/1.13 ivar
module load bowtie2
module load minimap2/2.22
module load edirect


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


export organisation=@fda.hhs.gov


echo Environment is set

