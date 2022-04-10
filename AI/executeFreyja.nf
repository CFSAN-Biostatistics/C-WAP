#! /usr/bin/env nextflow

// Input is grouped in subfolders to reduce the overhead due to thousands of slurm jobs taking 10s each.
inSubfolders = Channel.fromPath("$params.inpath/*", type: 'any', checkIfExists: true)


process freyja {
	time = '20 min'
	cpus = 1
	memory = '8 GB'
	errorStrategy = 'retry'
	maxRetries = 10
	
	publishDir "$params.outpath/", mode: 'copy', overwrite: true
	input:
		path('subfolder') from inSubfolders
	
	output:
		path('demix/*') into outfiles
		
	shell:
	"""
		unset PYTHONPATH
		freyjaBin=/projects/covidtrakr/software/miniconda3/envs/freyja-env/bin/freyja
		mkdir ./demix
		for fileIdx in \$(ls subfolder/depths/* | awk -F '/' '{print \$NF}' | awk -F '.' '{print \$1}'); do
			\$freyjaBin demix subfolder/variants/\$fileIdx.tsv subfolder/depths/\$fileIdx.tsv --output ./demix/\$fileIdx
		done
	"""
}

